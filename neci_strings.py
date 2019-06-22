def write_inp(fname, electrons, nwalkers, nmcyc=-1, time=None, write_wf=False, read_wf=False, write_rdm=False, shiftdamp=0.05, memoryfacspawn=20.0,
        read_rdm=False, granularity=1, cicj_promotion_hwhm_fraction=0, pops_core=0, write_core=False, read_core=False, two_rdm_only=False, stepsshift=5,
        ss_start=0, rdm_start=0, rdm_iters=0, rdm_energy_iters=None, walkcontgrow=False, seed=14, main_facs=(1,)*5, spawn_facs=(1,)*5, recv_facs=(1,)*5,
        hbrdm_offdiag_frac_occ_thresh=0.0, tau='', trial_start=0, pops_trial=0, niter_cicj_histogram=100, cicj_promotion_scale=1.0):
    if main_facs is not None: main_facs = ' '.join(map(str, main_facs))
    if spawn_facs is not None: spawn_facs = ' '.join(map(str, spawn_facs))
    if recv_facs is not None: recv_facs = ' '.join(map(str, recv_facs))
    with open(fname, 'w') as f:
        f.write('''title
system read noorder
electrons {}
symignoreenergies
freeformat
sym 0 0 0 0
nonuniformrandexcits 4ind-weighted
nobrillouintheorem
endsys
calc
methods
method vertex fcimc
endmethods
totalwalkers {}
memoryfacpart 10.0
memoryfacspawn {}
seed {}
startsinglepart 10
shiftdamp {}
truncinitiator
addtoinitiator 3.0
allrealcoeff
realspawncutoff 0.4
jump-shift
stepsshift {}
maxwalkerbloom 3
load-balance-blocks off
nmcyc {}
'''.format(electrons, nwalkers, memoryfacspawn, seed, shiftdamp, stepsshift, nmcyc))
        if write_rdm and rdm_iters>0: f.write('rdmsamplingiters {}\n'.format(rdm_iters))
        if read_wf or read_rdm: f.write('readpops\n')
        if len(tau): f.write('tau {}\n'.format(tau))
        if walkcontgrow: f.write('walkcontgrow\n')
        if time: f.write('time {}\n'.format(time))
        if read_core:
            if ss_start:
                f.write('semi-stochastic {}\n'.format(ss_start))
            else:
                f.write('semi-stochastic\n')
            f.write('read-core\n')
        elif pops_core > 0:
            if ss_start:
                f.write('semi-stochastic {}\n'.format(ss_start))
            else:
                f.write('semi-stochastic\n')
            f.write('pops-core {}\n'.format(pops_core))
        if pops_trial > 0:
            if trial_start:
                f.write('trial-wavefunction {}\n'.format(trial_start))
            else:
                f.write('trial-wavefunction\n')
            f.write('pops-trial {}\n'.format(pops_trial))
        f.write('''endcalc
integral
freeze 0 0
endint
logging
''')
        if write_core:
            assert pops_core > 0
            f.write('write-core\n')

        if write_rdm:
            if rdm_energy_iters is None: rdm_energy_iters = rdm_iters
            f.write('calcrdmonfly 3 {} {}\n'.format(rdm_start, rdm_energy_iters))
            if not two_rdm_only:
                f.write('nevpt2\n')
                f.write('hbrdm-ghost-granularity {}\n'.format(granularity))
                f.write('cicj-promotion-hwhm-fraction {}\n'.format(cicj_promotion_hwhm_fraction))
                f.write('cicj-promotion-scale {}\n'.format(cicj_promotion_scale))
                f.write('niter-cicj-histogram {}\n'.format(niter_cicj_histogram))
            f.write('write-spin-free-rdm\n')
            f.write('printonerdm\n')
            if main_facs is not None: f.write('rdm-main-size-fac {}\n'.format(main_facs))
            if spawn_facs is not None: f.write('rdm-spawn-size-fac {}\n'.format(spawn_facs))
            if recv_facs is not None: f.write('rdm-recv-size-fac {}\n'.format(recv_facs))
        if write_rdm or write_wf:
            f.write('binarypops\n')
        if read_rdm:
            f.write('readrdms\n')
        f.write('''endlog
end
''')
