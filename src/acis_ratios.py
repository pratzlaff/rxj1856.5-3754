import argparse
import astropy.io.fits
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pprint

from matplotlib import rc, rcParams
#rc('text', usetex=True)
#rcParams.update({'font.size': 14})
#rc('text', usetex=True)

import response
import util

def qdp_file(year, disp, args):
    return'{}/rxj1856_bbody_{}_{}_orders_5.qdp'.format(args.fitdir, year, disp)

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

def get_response_hrc(year, disp, args):
    arf = '{}/{}/{}_leg_{}1.arf'.format(args.hrcdir, disp, year, disp)
    rmf = '{}/{}/{}_leg_{}1.rmf'.format(args.hrcdir, disp, year, disp)

    bin_lo, bin_hi, resp = response.read_arf(arf)
    jnk1, jnk2, jnk3 = response.read_rmf(rmf)

    resp *= jnk3

    return 0.5*(bin_lo+bin_hi)[::-1], np.ma.masked_invalid(resp[::-1])

def get_response_acis(obsid):
    arf = './arfs/acis/{}_LEG_-1_garf.fits'.format(obsid)
    rmf = '/data/legs/rpete/flight/rmfs/ACIS-S-LEG_-1.rmf'

    bin_lo, bin_hi, resp = response.read_arf(arf)
    jnk1, jnk2, jnk3 = response.read_rmf(rmf)

    resp *= jnk3

    return 0.5*(bin_lo+bin_hi)[::-1], np.ma.masked_invalid(resp[::-1])

def main():
    parser = argparse.ArgumentParser(
        description='Compare ACIS-S spectra with 2001 HRC-S combined spectrum.'
    )
    parser.add_argument('-y', '--ylabel', help='Y axis label', default=r'\textrm{New / Old}')
    parser.add_argument('-t', '--title', help='Plot title', default=r'\textrm{Rate Ratio}')
    parser.add_argument('-o', '--outfile', help='Output file name')
    parser.add_argument('--hrcdir', default='./spectra/qe_N0015_qeu_N0013', help='Directory containing HRC-S spectra.')
    parser.add_argument('--fitdir', default='./fits/qe_N0015_qeu_N0013', help='Directory containing HRC-S fits.')
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')

    args = parser.parse_args()

    min, max, inc = 20., 50, 10.
    n = int((max-min)/inc)
    w1 = np.arange(n+1)*inc+min
    w2 = w1 + inc

    d = {}
    for year in 2001,:
        d[year] = {}
        for disp in 'm':
            wav, src, hdr = read_pha('{}/{}/{}_leg_{}1.pha'.format(args.hrcdir, disp, year, disp))
            bg, backscal = read_bkg('{}/{}/{}_leg_{}1_bkg.pha'.format(args.hrcdir, disp, year, disp, disp))
            resp_wav, resp = get_response_hrc(year, disp, args)
            frac_wav, frac = qdp_ratios(qdp_file(year, disp, args), args)

            resp = np.interp(wav, resp_wav, resp)
            frac = np.interp(wav, frac_wav, frac)

            d[year] = {
                'wav' : wav,
                'src' : src,
                'hdr' : hdr,
                'bg' : bg,
                'backscal' : backscal,
                'response' : resp,
                'frac' : frac,
            }

    obsids = [ 13198, 14267 ]
    row = 0 # TG_M = -1
    for obsid in obsids:
        pha2 = glob.glob('./data/{}/tg_reprocess/*_pha2.fits'.format(obsid))[0]
        data, hdr = util.read_pha2(pha2)
        wav = 0.5 * (data['bin_lo'][row]+data['bin_hi'][row])[::-1]
        resp_wav, resp = get_response_acis(obsid)
        resp = np.interp(wav, resp_wav, resp)

        d[obsid] = {
            'wav' : wav,
            'src' : data['counts'][row][::-1],
            'bg' : (data['background_up'][row] + data['background_down'][row])[::-1],
            'backscal' : np.zeros_like(wav) + hdr['backscup'] + hdr['backscdn'],
            'response' : resp,
            'frac' : np.ones_like(wav),
            'date-obs' : hdr['date-obs'][:10],
        }
            
    for key in d:
        wav = d[key]['wav']
        src = d[key]['src']
        bg = d[key]['bg']
        backscal = d[key]['backscal']
        resp = d[key]['response']

        for quant in 'rate', 'rate_err', 'flux', 'flux_err':
            d[key][quant] = np.zeros_like(w1)

        for i in range(d[key]['rate'].size):
            ii = np.where((wav>=w1[i]) & (wav<w2[i]))

            r, rerr, f, ferr = flux_summed(wav[ii], src[ii], bg[ii], backscal[ii], resp[ii], hdr, fraction=frac[ii])

            d[key]['rate'][i] = r
            d[key]['rate_err'][i] = rerr
            d[key]['flux'][i] = f
            d[key]['flux_err'][i] = ferr

    plot_dims = (4, 1)
    figsize = (8.5, 11)
    if args.pdf:
        #rcParams.update({'font.size': 14})
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = figsize)

    wav = 0.5*(w1+w2)

    for i in range(len(obsids)):
        obsid = obsids[i]
        row = i % plot_dims[0]
        plt.subplot2grid(plot_dims, (row,0))
        plt.title('RX J1856.5-3754: ACIS-S / HRC-S (2001) - {}: {}'.format(obsid, d[obsid]['date-obs']), loc='right')
        plt.ylabel('Flux Ratio')
        

        x = wav
        r1 = d[obsid]['flux']
        r1err = d[obsid]['flux_err']
        r2 = d[13198]['flux']
        r2err = d[13198]['flux_err']

        r = r1/r2
        rerr = r * np.sqrt( (r1err/r1)**2 + (r2err/r2)**2 )

        plt.errorbar(x, r, yerr=rerr)

        if i == len(obsids)-1 or (i % plot_dims[0])==plot_dims[0]-1:
            plt.xlabel(r'$\lambda\;(\AA)$')

            plt.tight_layout()
            if args.pdf:
                pdf.savefig(fig)
            else:
                plt.show()
            plt.clf()

    if args.pdf:
        pdf.close()

    plt.close()

if __name__ == '__main__':
    main()
