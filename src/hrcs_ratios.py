import astropy.io.fits
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
#rc('text', usetex=True)
#rcParams.update({'font.size': 14})
#rc('text', usetex=True)

ANGSTROM, LAMBDA = "Åλ"

import response
import qdp

def qdp_file(year, disp, args):
    return f'{args.fitdir}/rxj1856_bbody_{year}_{disp}_orders_5.qdp'

def read_qdp(qdpfile, args):
    data = np.genfromtxt(qdpfile, skip_header=3, loose=True, unpack=True)
    n = int((data.shape[1]-1)/2)
    return data[0,:n], data[5:,:n]

def qdp_ratios(qdpfile, args):
    wav, orders = read_qdp(qdpfile, args)
    ratio = orders[0] / orders.sum(axis=0)
    return wav[::-1], ratio[::-1]

def flux_summed(wav, src, bg, backscal, resp, hdr, fraction=1):
    net = src - bg / backscal
    net_var = src + bg / backscal / backscal

    rate = net / hdr['exposure']
    rate_err = np.sqrt(net_var) / hdr['exposure']

    rate *= fraction
    rate_err *= fraction

    energy = 12.39854 / wav
    mult = energy * 1.602e-9 / resp
    flux = rate * mult
    flux_err = rate_err * mult

    rate = rate.sum()
    rate_err = np.sqrt((rate_err**2).sum())

    flux = flux.sum()
    flux_err = np.sqrt((flux_err**2).sum())

    return rate, rate_err, flux, flux_err

def read_pha(phafile):
    hdulist = astropy.io.fits.open(phafile)
    data = hdulist[1].data
    bin_lo, bin_hi, src = data.field('bin_lo'), data.field('bin_hi'), data.field('counts')
    hdr = hdulist[1].header
    hdulist.close()

    return 0.5*(bin_lo+bin_hi)[::-1], src[::-1], hdr

def read_bkg(bgfile):
    hdulist = astropy.io.fits.open(bgfile)
    data = hdulist[1].data
    bg, backscal = data.field('counts'), data.field('backscal')
    hdulist.close()

    return bg[::-1], backscal[::-1]

def get_response(year, disp, args):
    arf = f'{args.hrcdir}/{disp}/{year}_leg_{disp}1.arf'
    rmf = f'{args.hrcdir}/{disp}/{year}_leg_{disp}1.rmf'

    bin_lo, bin_hi, resp = response.read_arf(arf)
    jnk1, jnk2, jnk3 = response.read_rmf(rmf)

    resp *= jnk3

    return 0.5*(bin_lo+bin_hi)[::-1], np.ma.masked_invalid(resp[::-1])

def main():
    parser = argparse.ArgumentParser(
        description='Compare (Type I) spectra, plotting new/old.'
    )
    parser.add_argument('-y', '--ylabel', help='Y axis label', default=r'\textrm{New / Old}')
    parser.add_argument('-t', '--title', help='Plot title', default=r'\textrm{Rate Ratio}')
    parser.add_argument('-o', '--outfile', help='Output file name')
    parser.add_argument('--hrcdir', default='./spectra/qeu_N0014', help='Directory containing HRC-S spectra.')
    parser.add_argument('--fitdir', default='./fits/qeu_N0014', help='Directory containing HRC-S fits.')
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')

    args = parser.parse_args()

    min, max, inc = 20., { 'm' : 82, 'p' : 82 }, 2.
    n, w1, w2 = {}, {}, {}
    for o in max:
        n[o] = int((max[o]-min)/inc)
        w1[o] = np.arange(n[o]+1)*inc+min
        w2[o] = w1[o] + inc

    d = {}
    for year in 2001, 2013, 2019:
        d[year] = {}
        for disp in 'm', 'p':
            wav, src, hdr = read_pha(f'{args.hrcdir}/{disp}/{year}_leg_{disp}1.pha')
            bg, backscal = read_bkg(f'{args.hrcdir}/{disp}/{year}_leg_{disp}1_bkg.pha')
            resp_wav, resp = get_response(year, disp, args)
            frac_wav, frac = qdp_ratios(qdp_file(year, disp, args), args)

            resp = np.interp(wav, resp_wav, resp)
            frac = np.interp(wav, frac_wav, frac)

            d[year][disp] = {
                'wav' : wav,
                'src' : src,
                'hdr' : hdr,
                'bg' : bg,
                'backscal' : backscal,
                'response' : resp,
                'frac' : frac,
            }

            for key in 'rate', 'rate_err', 'flux', 'flux_err':
                d[year][disp][key] = np.zeros_like(w1[disp])

            for i in range(d[year][disp]['rate'].size):
                ii = np.where((wav>=w1[disp][i]) & (wav<w2[disp][i]))

                r, rerr, f, ferr = flux_summed(wav[ii], src[ii], bg[ii], backscal[ii], resp[ii], hdr, fraction=frac[ii])

                d[year][disp]['rate'][i] = r
                d[year][disp]['rate_err'][i] = rerr
                d[year][disp]['flux'][i] = f
                d[year][disp]['flux_err'][i] = ferr

    formats = { 'm' : 'ro--', 'p' : 'k>-' }
    labels = { 'm' : 'TG_M=-1', 'p' : 'TG_M=+1' }
    offsets = { 'm' : -0.1, 'p' : 0.1 }
    colors = { 'm' : 'r', 'p' : 'k' }

    ylims = {
        2013 : { 'rate' : [0.5, 1.2], 'flux' : [0.6, 1.23] },
        2019 : { 'rate' : [0.55, 1.0], 'flux' : [0.77, 1.24] }
    }

    plot_dims = (2, 1)
    figsize = (8.5, 11)
    if args.pdf:
        rcParams.update({'font.size': 14})
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = figsize)

    wav = {}
    for o in max:
        wav[o] = 0.5*(w1[o]+w2[o])
    for year in 2013, 2019:

        plt.subplot2grid(plot_dims, (0,0))
        plt.title(f'RX J1856.5-3754: {year} / 2001')
        plt.ylabel('Rate Ratio')
        global LAMBDA, ANGSTROM
        plt.xlabel(f'{LAMBDA} {ANGSTROM})')
        x = {}
        for disp in d[year]:
            x[disp] = wav[disp] + offsets[disp]
            
            r1 = d[year][disp]['rate']
            r1err = d[year][disp]['rate_err']
            r2 = d[2001][disp]['rate']
            r2err = d[2001][disp]['rate_err']

            r = r1/r2
            rerr = r * np.sqrt( (r1err/r1)**2 + (r2err/r2)**2 )

            plt.errorbar(x[disp], r, yerr=rerr, fmt=formats[disp], label=labels[disp])
        plt.ylim(ylims[year]['rate'])
        plt.legend()
                
        plt.subplot2grid(plot_dims, (1,0))
        plt.ylabel('Flux Ratio')
        plt.xlabel('Wavelength')
        plt.xlabel(f'{LAMBDA} {ANGSTROM})')
        for disp in d[year]:
            f1 = d[year][disp]['flux']
            f1err = d[year][disp]['flux_err']
            f2 = d[2001][disp]['flux']
            f2err = d[2001][disp]['flux_err']

            r = f1/f2
            rerr = r * np.sqrt( (f1err/f1)**2 + (f2err/f2)**2 )

            plt.errorbar(x[disp], r, yerr=rerr, fmt=formats[disp], label=labels[disp])
        plt.ylim(ylims[year]['flux'])
        plt.legend()

        plt.tight_layout()
        if args.pdf:
            pdf.savefig(fig)
        else:
            plt.show()
        plt.clf()
        

    # plt.plot(0.5 * (w1+w2), ratio)
    # plt.tight_layout()
    # plt.ylim(0.6, 1)
    # plt.show()
    # e1, qe1 = read_qe(args.old, 1)
    # e2, qe2 = read_qe(args.new, 1)

    # ratio = qe2 / np.interp(e2, e1, qe1)
    # wav = 12.39854 / e2
    # plt.plot(wav, ratio)
    # plt.xlabel(r'\textrm{Wavelength (\AA)}')
    # plt.ylabel(args.ylabel)
    # plt.title(args.title)

    # plt.tight_layout()

    # if (args.outfile): plt.savefig(args.outfile)
    # else: plt.show()

    # plt.close()

    if args.pdf:
        pdf.close()

if __name__ == '__main__':
    main()
