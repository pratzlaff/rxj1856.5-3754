import glob

from sherpa.astro.ui import *

obsid=13198

order=-1; row=1
#order=1; row=4

pha2=glob.glob('/data/legs/rpete/flight/rxj1856.5-3754/data/{}/tg_reprocess/*_pha2.fits'.format(obsid))[0]

set_default_id(row)
load_pha(pha2, use_errors=True)

group_counts(1)

set_stat('wstat')

garf='/data/legs/rpete/flight/rxj1856.5-3754/arfs/acis/{}_LEG_{}_garf.fits'.format(obsid, order)
rmf='/data/legs/rpete/flight/rmfs/ACIS-S-LEG_{}.rmf'.format(order)
load_arf(garf)
load_rmf(rmf)

wlo, whi = 10, 52

set_analysis('wave')
ignore()
notice(wlo, whi)

set_xsabund("wilm")
set_source(xstbabs.tbabs * xsbbody.bb)
tbabs.nh = 0.01
freeze(tbabs.nh)
bb.kt = 0.065
bb.norm = 2e-4

fit()

set_ylog()
ungroup()
group_snr(5)
plot_fit_ratio()
