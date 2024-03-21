import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scienceplots


plt.style.use(['science', 'notebook'])

hbarc = 197.3269804
Gamma0 = 1
nfil=6


Path=input('Enter the relative or absolute path to directory where the "*.res" and for005 (for the dynamic Run) files are located: \n')
# print(Path)

with open(Path+'/for005','r') as f:
    lines=f.readlines()
    word1='L_val'
    word2='M_val'
    for line in lines:
        # check if string present on a current line
        if line.find('force') != 1:
            data = line.split(',')
            for d in data:
                if d.find('force') != -1:
                    force = d.split('=')[1]

        if line.find(word1) != -1:
            # print(word, 'string exists in file')
            # print('Line Number:', lines.index(line))
            # print('Line:', line.split(','))
            data = line.split(',')
            for d in data:
                if d.find(word1) !=-1:
                    # print(d)
                    L_val = int(d.split('=')[1])
                elif d.find(word2) != -1:
                    # print(d)
                    M_val = int(d.split('=')[1])
            exit

print(f'L = {L_val}, M = {M_val}, force = {force}')


def FileNameFromL(L):
    filenames={0:'monopoles.res',1:'dipoles.res',2:'quadrupoles.res',3:'octupoles.res',4:'hexadecapoles.res',5:'diatriacontapoles.res'}
    return filenames[L]


if Path[-1] == '/':
    OutPutFile = 'Strength_'+FileNameFromL(L_val).replace('res', 'out')
    ResFile = Path+FileNameFromL(L_val)
    ExtFieldFile=Path+'extfield.res'
else:
    OutPutFile = '/Strength_'+FileNameFromL(L_val).replace('res', 'out')
    ResFile = Path+'/'+FileNameFromL(L_val)
    ExtFieldFile=Path+'/extfield.res'
print(ResFile,ExtFieldFile)
print(OutPutFile)


if L_val == 0:
    ylabel = r'S(E) (fm$^4$/MeV)'
    legend='Monopole Boost'
elif L_val == 1:
    ylabel = r'S(E) (fm$^6$/MeV)'
    legend='Dipole Boost'
elif L_val == 2:
    ylabel = r'S(E) (fm$^4$/MeV)'
    legend='Quadrupole Boost'
elif L_val == 3:
    ylabel = r'S(E) (fm$^6$/MeV)'
    legend='Octupole Boost'
elif L_val == 4:
    ylabel = r'S(E) (fm$^8$/MeV)'
    legend='Hexadecapole Boost'
elif L_val == 5:
    ylabel = r'S(E) (fm$^{{{10}}}$/MeV)'
    legend='Diatriaconta (32) Boost'



def Read(resFile, extFile):

    with open(resFile, 'r') as f:
        lines = f.readlines()[1:]
        time1 = np.zeros((len(lines)), dtype=float)
        Signal1 = np.zeros((len(lines)), dtype=float)
        for i, line in enumerate(lines):
            data = line.split()
            time1[i] = float(data[0])
            Signal1[i] = float(data[1])
        amplitude = float(lines[0].split()[2])
        Signal1=Signal1-Signal1[0]
        
    with open(extFile, 'r') as f:
        lines = f.readlines()[1:]
        time2 = np.zeros((len(lines)), dtype=float)
        Signal2 = np.zeros((len(lines)), dtype=float)
        for i, line in enumerate(lines):
            data = line.split()
            time2[i] = float(data[0])
            Signal2[i] = float(data[1])
        amplitude = float(lines[0].split()[2])
        L = int(lines[0].split()[3])
        M = int(lines[0].split()[4])
        Signal2=Signal2-Signal2[0]
    # print(time1,Signal1)
    return time1, Signal1, time2, Signal2,L,M,amplitude


def Filter1(Gamma0, t, hbc):
    return np.exp((-Gamma0*t)/(2.0*hbc))

def Filter2(nfil,t):
    return (np.cos(0.5*np.pi*t/t[-1]))**nfil

# plt.plot(time,Quad_amp*Filter(Gamma0,time,hbarc))
# plt.show()

def Fourier(y, t):
    Y = np.fft.rfft(y, norm='backward')*2
    freq = np.fft.rfftfreq(len(t), t[1]-t[0])
    # print(t[1]-t[0], freq*hbarc*2*np.pi)
    return freq*hbarc*2*np.pi, Y, Y.imag



time1, Moment1, time2, Moment2, L, M, ampl = Read(ResFile, ExtFieldFile)

division_const = ampl*hbarc*np.pi

Moment1_Exp_Filter = Moment1*Filter1(Gamma0,time1,hbarc) 
Moment1_Cos_Filter = Moment1*Filter2(nfil,time1)

E1_Exp_fil, Spectrum_c1_exp_fil, Spectrum1_exp_fil = Fourier(Moment1_Exp_Filter, time1)
E1_cos_fil, Spectrum_c1_cos_fil, Spectrum1_cos_fil = Fourier(Moment1_Cos_Filter, time1)

Spectrum1_cos_fil = Spectrum1_cos_fil/division_const
Spectrum1_exp_fil = Spectrum1_exp_fil/division_const
Spectrum_c1_cos_fil = Spectrum_c1_cos_fil/division_const
Spectrum_c1_exp_fil = Spectrum_c1_exp_fil/division_const


Moment2_Exp_Filter = Moment2*Filter1(Gamma0,time2,hbarc)
Moment2_Cos_Filter = Moment2*Filter2(nfil,time2)

E2_Exp_fil, Spectrum_c2_exp_fil, Spectrum2_exp_fil = Fourier(Moment2_Exp_Filter, time2)
E2_cos_fil, Spectrum_c2_cos_fil, Spectrum2_cos_fil = Fourier(Moment2_Cos_Filter, time2)

Spectrum2_cos_fil = Spectrum2_cos_fil/division_const
Spectrum2_exp_fil = Spectrum2_exp_fil/division_const
Spectrum_c2_cos_fil = Spectrum_c2_cos_fil/division_const
Spectrum_c2_exp_fil = Spectrum_c2_exp_fil/division_const



with open(Path+OutPutFile,'w') as f:
    f.write('{:<10} \t {:<10} \t {:<10}  \n'.format('Energy(MeV)', 'Strength (Exp Filter)', 'Strength (Cos Filter)  Units = ['+ylabel[4:]+']'))
    for E,S1,S2 in zip(E1_Exp_fil,Spectrum1_exp_fil,Spectrum1_cos_fil):
        f.write('{:<12.2f} \t {:<18.7f} \t {:<12.7f} \n'.format(E,S1,S2))

print('\n -------------------------------------------------')
print(
    f'Strength function calculated and written in the file = {Path+OutPutFile}')
print('--------------------------------------------------')



