import bisect, os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import matplotlib
'''
new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)
'''
def get_loop_loci(loops):

    loci = {}
    for l in loops:
        c1, s1, e1, c2, s2, e2 = l
        if not c1 in loci:
            loci[c1] = set()
        if not c2 in loci:
            loci[c2] = set()
        loci[c1].add((s1, e1))
        loci[c2].add((s2, e2))

    return loci

def parseLoops(loopfil):

    loops = []
    with open(loopfil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            c1 = parse[0].lstrip('chr')
            c2 = parse[3].lstrip('chr')
            loops.append([c1, int(parse[1]), int(parse[2]), c2, int(parse[4]), int(parse[5])])
    
    return loops

def parseChIP(fil):
    
    data = {}
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0].lstrip('chr')
            if '_' in chrom:
                continue
            if chrom in ['Y','M']:
                continue
            if chrom in data:
                data[chrom].append([int(parse[1]), int(parse[2])])
            else:
                data[chrom] = [[int(parse[1]), int(parse[2])]]
    for chrom in data:
        data[chrom].sort()
    
    return data

def centerLoci(loci, chip, halfr=50000, bins=50):

    arrays = []
    for c in loci:
        if not c in chip:
            continue
        bychrom = chip[c]
        for p in loci[c]:
            center = (p[0] + p[1]) // 2
            r = [c, center-halfr, center+halfr]
            data = np.zeros(r[2]-r[1])
            idx = max(0, bisect.bisect(bychrom, r[1:])-1)
            for q in bychrom[idx:]:
                if q[1] <= r[1]:
                    continue
                if q[0] >= r[2]:
                    break
                s = q[0]-r[1]
                if s < 0:
                    s = 0
                e = q[1]-r[1]
                data[s:e] += 1
            idx = np.linspace(0,r[2]-r[1],bins+1).astype(int)
            arr = []
            for i, j in zip(idx[:-1], idx[1:]):
                arr.append(data[i:j].mean())
            arr = np.r_[arr]
            '''
            if arr.sum() > 0:
                arr = arr / arr.mean()
            '''
            arrays.append(arr)
    
    arr = np.r_[arrays].mean(axis=0)
    num = len(arrays)
    
    return arr, num

loops1 = parseLoops(sys.argv[1]) # bwa unique, bedpe
loop_loci1 = get_loop_loci(loops1)
loops2 = parseLoops(sys.argv[2]) # chromap unique, bedpe
loop_loci2 = get_loop_loci(loops2)
loops3 = parseLoops(sys.argv[3]) # common, bedpe
loop_loci3 = get_loop_loci(loops3)
chip = parseChIP(sys.argv[4]) # CTCF ChIP-Seq peaks, bed
cell = sys.argv[1].split('.')[0]

halfr = 300000
arr1, num1 = centerLoci(loop_loci1, chip, halfr=halfr, bins=51)
arr2, num2 = centerLoci(loop_loci2, chip, halfr=halfr, bins=51)
arr3, num3 = centerLoci(loop_loci3, chip, halfr=halfr, bins=51)
fig = plt.figure(figsize=(2.2, 1.8))
plt.rc('font', size=6)
ax = fig.add_subplot(111)

l3, = ax.plot(arr3, color='#FDC086', lw=1.5)
l1, = ax.plot(arr1, color='#7FC97F', lw=1.5)
l2, = ax.plot(arr2, color='#BEAED4', lw=1.5)
ax.set_xticks([0,25,50])
ax.set_xticklabels(['-300K','loop anchor','+300K'], fontsize=6) 
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel('Normalized number of CTCF peaks', fontsize=6)

ax.legend([l3, l1, l2], ['common', 'bwa unique', 'chromap unique'], frameon=False, fontsize=6)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)

ax.xaxis.set_tick_params(width=1, labelsize=6, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=6, pad=2)
ax.set_axisbelow(True)

plt.savefig('{0}.chip-enrich.svg'.format(cell), dpi=300, bbox_inches='tight')
plt.close()
