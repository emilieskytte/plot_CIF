import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
import re
import os

#filenames = [i.split("\\")[-1] for i in sys.argv[1:-1]]
filenames = sys.argv[1:-1]

no_atoms = int(sys.argv[-1])

path =  os.path.dirname(sys.argv[1])

def get_ADPs(file, no_atoms, *args, **kwargs):
    with open(file) as ofile:
        rfile = ofile.readlines()
    
    start = rfile.index('  _atom_site_aniso_U_12\n')
    end = start+1+no_atoms
    
    ip = rfile[start+1:end]

    ADPs = {}
    for a in ip:
        atoms = a.split()[0]
        adp_picker = re.findall('(\d+.\d+)\((\d+)\)', a)[0:3]
        ADPs[atoms] =  [(float(adp_picker[i][0]), float(adp_picker[i][1])*10**-len(adp_picker[i][0].split('.')[1])) for i in range(0,3)]
        
    
    return ADPs
    
    
def get_occupancy(file, no_atoms, *args, **kwargs):
    with open(file) as ofile:
        rfile = ofile.readlines()
    
    start = rfile.index('  _atom_site_disorder_group\n')
    end = start+1+no_atoms
    
    ip = rfile[start+1:end]
    atom_str = [i.split() for i in ip]
    
    atom = [i[0] for i in atom_str]
    occ_str = [i[7] for i in atom_str]

    occ = {}
    for a in zip(atom, occ_str):
        occ_float = re.findall('(\d+.\d+)(\d+)', a[1])
      
        if len(occ_float) == 0:
            continue  
        occ[a[0]] = [float(occ_float[0][0]), float(occ_float[0][1])*10**-len(occ_float[0][0].split('.')[1])]

    return occ
    
    

def get_cell(file, *args, **kwargs):
    with open(file) as ofile:
        rfile = ofile.readlines()
    
    length, angle, volume = [], [], []
    
    for i in rfile:
        length += re.findall('_cell_length_(.{1})\s*(\d+.\d+).(\d+)', i)
        angle += re.findall('_cell_angle_(.+)\s*(\d+.\d+).(\d+)', i)
        volume += re.findall('_cell_volume\s*(\d+.\d+).(\d+)', i)

    cell = {}
    for a in [*length,*angle,('volume', *volume[0])]:
        if a[0] == 'b':
            continue
        elif a[0] == 'c':
            continue
        cell[a[0]] = [float(a[1]), float(a[2])*10**-len(a[1].split('.')[1])]
    
    return cell
    
    
def all_files(file_list, func,  *args, **kwargs):
    rv = func(file_list[0],no_atoms, *args,  **kwargs)
    for key, value in rv.items():
        value = [[v] for v in value]
        rv[key] = value
        for f in file_list[1:]:
            new_element = func(f, no_atoms, *args, **kwargs)
            for i, ne in enumerate(rv[key]):
                rv[key][i].append(new_element[key][i])
            
    return rv

def get_temp(file_list):
    temp = []
    for f in file_list:
        t = re.findall('.*_(\d+)K_.*', f)
        temp.append(float(t[0]))
    return temp

def get_Ueq(adp_dict):
    ueq_dict = {}
    for key, value in adp_dict:
            ueq_dict[key] = [np.mean(zip(value)), np.std(zip(value))]
    return ueq_dict

#print(get_ADPs(filenames[1], no_atoms))
#print(all_files(filenames, get_occupancy, no_atoms))
T = get_temp(filenames)


cmap = [cm.plasma(x) for x in np.linspace(0.0, 1.0, no_atoms)]

### plot ADPs
fig, axs = plt.subplots(no_atoms)
fig.set_size_inches(7,4*no_atoms)
i=0
for atom, value in all_files(filenames, get_ADPs, no_atoms).items():
    axp = axs[i]

    zip_adp = [list(zip(*uij)) for uij in value]

    #m,b,merr,berr = lin_reg(T[1::],np.array(u[e][1::]))
    #D, Derr = debyeT(m,merr,mass[0])
    ueq = [np.mean(v) for v in list(zip(*[i[0] for i in zip_adp]))]
    ueq_err =[np.mean(v) for v in list(zip(*[i[1] for i in zip_adp]))]

    axp.errorbar(T, ueq, yerr=ueq_err, fmt='o-',c=cmap[i], mec='k', capsize=2, label=atom)
    
    #axp.plot(x, m*x+b,':',c=cmap[0], label=r'$\theta_D$ = {:.0f}({:.0f})'.format(D, Derr*10**0))
    ulabel = [r'$U_{11}$',r'$U_{22}$',r'$U_{33}$']
    for ii, uij in enumerate(zip_adp):
        axp.errorbar(T, uij[0], yerr=uij[1], fmt='o-', c=cmap[i], mec='k', alpha = 0.2, capsize=2, label=ulabel[ii])
    

    axp.legend()
    axp.set_xlabel('Temperature [K]')  
    axp.set_ylabel('$U_{eq}$ [Å$^2$]')  

    i+=1
fig.savefig(path + r'\fig_adp.png' ,dpi=200,bbox_inches = "tight")




### plot occupancy
no_occ = len(get_occupancy(filenames[0], no_atoms))
fig, axs = plt.subplots(no_occ)
fig.set_size_inches(7,4*no_occ)
i=0
for atom, value in all_files(filenames, get_occupancy, no_atoms).items():
    axp = axs
    if no_occ != 1:
        axp = axs[i]

    axp.errorbar(T, value[0], yerr=value[1],fmt='o-',c=cmap[i], mec='k', capsize=2, label=atom)
    
    axp.legend()
    axp.set_xlabel('Temperature [K]')  
    axp.set_ylabel('Occupancy [frac.]')  
    i+=1
fig.savefig(path + r'\fig_occ.png',dpi=200,bbox_inches = "tight")





### plot cell
no_cell = len(get_cell(filenames[0]))
fig, axs = plt.subplots(no_cell)
fig.set_size_inches(7,4*no_cell)
i=0
for key, value in all_files(filenames, get_cell).items():
    axp = axs[i]

    
    
    p,cov=np.polyfit(T,(np.array(value[0])-value[0][0])/value[0][0]*100,1,cov=True)
    m, b = p
    merr = np.sqrt(cov[0][0])
    berr = np.sqrt(cov[1][1])

    axp.errorbar(T, value[0], yerr=value[1], fmt='o-',c=cmap[i], mec='k', capsize=2,  label=r'$\Delta${}% = {:.4f}({:.0f}) $\times$ T - {:.2f}({:.0f})'.format(key, m, merr*10**4, abs(b), berr*10**2))
    
    axp.legend()
    axp.set_xlabel('Temperature [K]')  
    axp.set_ylabel('Unit cell parameter [Å]')  
    if len(key) > 1:
        axp.set_ylabel('Unit cell angle [$^o$]') 
    if key == 'volume':
        axp.set_ylabel('Unit cell volume [Å$^3$]') 
    i+=1
#savename = os.path.join(path,'fig_cell.png')
fig.savefig(path + r'\fig_cell.png',dpi=200,bbox_inches = "tight")

















