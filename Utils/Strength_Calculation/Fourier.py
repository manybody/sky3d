import numpy as np



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
        if line.find('force') != -1:
            data = line.split(',')
            for d in data:
                if d.find('force') != -1:
                    force = d.split('=')[1]
        
        if line.find('ampl_ext') != -1:
            data = line.split(',')
            print(line)
            for d in data:
                if d.find('ampl_ext') != -1:
                    ampl = float(d.split('=')[1].replace('D','e'))

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

print(f'L = {L_val}, M = {M_val}, force = {force}, Ext. Amplitude = {ampl:f}')
def FileNameFromL(L):
    filenames={0:'monopoles.res',1:'dipoles.res',2:'quadrupoles.res',3:'octupoles.res',4:'hexadecapoles.res',5:'diatriacontapoles.res'}
    return filenames[L]


Coupled = input('Do you want to analyse a different multipole file to that of applied boost? \n(Leave Empty of write False if answer is No)')
if Coupled:
    L_val = int('Which Multipole File ? Enter a Integer value to corresponding L (0:monopoles, 1:dipoles .. so on)')
    print(f'File Chosen = {FileNameFromL(L_val)}')

if Path[-1] == '/':
    OutPutFile_IS = 'IS_Strength_'+FileNameFromL(L_val).replace('res', 'out')
    OutPutFile_IV = 'IV_Strength_'+FileNameFromL(L_val).replace('res', 'out')
    ResFile = Path+FileNameFromL(L_val)
    ExtFieldFile = Path+'extfield.res'
else:
    OutPutFile_IS = '/IS_Strength_'+FileNameFromL(L_val).replace('res', 'out')
    OutPutFile_IV = '/IV_Strength_'+FileNameFromL(L_val).replace('res', 'out')
    ResFile = Path+'/'+FileNameFromL(L_val)
    ExtFieldFile = Path+'/extfield.res'
print(ResFile, ExtFieldFile)
print(OutPutFile_IS, OutPutFile_IV)


def GetUnits(L_val, kind):
    if L_val == 0:
        ylabel = r'S(E) (fm$^4$/MeV)'
        legend = 'Monopole Boost'
    elif L_val == 1:
        if kind == 'IS':
            ylabel = r'S(E) (fm$^6$/MeV)'
        elif kind == 'IV':
            ylabel = r'S(E) (fm$^2$/MeV)'
        legend = 'Dipole Boost'
    elif L_val == 2:
        ylabel = r'S(E) (fm$^4$/MeV)'
        legend = 'Quadrupole Boost'
    elif L_val == 3:
        ylabel = r'S(E) (fm$^6$/MeV)'
        legend = 'Octupole Boost'
    elif L_val == 4:
        ylabel = r'S(E) (fm$^8$/MeV)'
        legend = 'Hexadecapole Boost'
    elif L_val == 5:
        ylabel = r'S(E) (fm$^{{{10}}}$/MeV)'
        legend = 'Diatriaconta (32) Boost'
    return ylabel, legend


def Read(resFile, extFile):

    with open(resFile, 'r') as f:
        lines = f.readlines()[1:]
        time1 = np.zeros((len(lines)), dtype=float)
        Signal_IS = np.zeros((len(lines)), dtype=float)
        Signal_IV = np.zeros((len(lines)), dtype=float)
        for i, line in enumerate(lines):
            data = line.split()
            time1[i] = float(data[0])
            Signal_IS[i] = float(data[1])
            Signal_IV[i] = float(data[2])
        Signal_IS = Signal_IS-Signal_IS[0]
        Signal_IV = Signal_IV-Signal_IV[0]

    with open(extFile, 'r') as f:
        lines = f.readlines()[1:]
        time2 = np.zeros((len(lines)), dtype=float)
        Signal2 = np.zeros((len(lines)), dtype=float)
        for i, line in enumerate(lines):
            data = line.split()
            time2[i] = float(data[0])
            Signal2[i] = float(data[1])
        Signal2 = Signal2-Signal2[0]
    # print(time1,Signal1)
    return time1, Signal_IS, Signal_IV, time2, Signal2


def Filter1(Gamma0, t, hbc):
    return np.exp((-Gamma0*t)/(2.0*hbc))


def Filter2(nfil, t):
    return (np.cos(0.5*np.pi*t/t[-1]))**nfil

# plt.plot(time,Quad_amp*Filter(Gamma0,time,hbarc))
# plt.show()


def Fourier(y, t):
    Y = np.fft.rfft(y, norm='backward')*2
    freq = np.fft.rfftfreq(len(t), t[1]-t[0])
    # print(t[1]-t[0], freq*hbarc*2*np.pi)
    return freq*hbarc*2*np.pi, Y, Y.imag


time1, Moment1_IS, Moment1_IV, time2, Moment2 = Read(ResFile, ExtFieldFile)

division_const = ampl*hbarc*np.pi

Moment1_IS_Exp_Filter = Moment1_IS*Filter1(Gamma0, time1, hbarc)
Moment1_IS_Cos_Filter = Moment1_IS*Filter2(nfil, time1)
Moment1_IV_Exp_Filter = Moment1_IV*Filter1(Gamma0, time1, hbarc)
Moment1_IV_Cos_Filter = Moment1_IV*Filter2(nfil, time1)

E1_Exp_fil_IS, Spectrum_c1_exp_fil_IS, Spectrum1_exp_fil_IS = Fourier(
    Moment1_IS_Exp_Filter, time1)
E1_cos_fil_IS, Spectrum_c1_cos_fil_IS, Spectrum1_cos_fil_IS = Fourier(
    Moment1_IS_Cos_Filter, time1)
E1_Exp_fil_IV, Spectrum_c1_exp_fil_IV, Spectrum1_exp_fil_IV = Fourier(
    Moment1_IV_Exp_Filter, time1)
E1_cos_fil_IV, Spectrum_c1_cos_fil_IV, Spectrum1_cos_fil_IV = Fourier(
    Moment1_IV_Cos_Filter, time1)


Spectrum1_cos_fil_IS = Spectrum1_cos_fil_IS/division_const
Spectrum1_exp_fil_IS = Spectrum1_exp_fil_IS/division_const
Spectrum_c1_cos_fil_IS = Spectrum_c1_cos_fil_IS/division_const
Spectrum_c1_exp_fil_IS = Spectrum_c1_exp_fil_IS/division_const

Spectrum1_cos_fil_IV = Spectrum1_cos_fil_IV/division_const
Spectrum1_exp_fil_IV = Spectrum1_exp_fil_IV/division_const
Spectrum_c1_cos_fil_IV = Spectrum_c1_cos_fil_IV/division_const
Spectrum_c1_exp_fil_IV = Spectrum_c1_exp_fil_IV/division_const


Moment2_Exp_Filter = Moment2*Filter1(Gamma0, time2, hbarc)
Moment2_Cos_Filter = Moment2*Filter2(nfil, time2)

E2_Exp_fil, Spectrum_c2_exp_fil, Spectrum2_exp_fil = Fourier(
    Moment2_Exp_Filter, time2)
E2_cos_fil, Spectrum_c2_cos_fil, Spectrum2_cos_fil = Fourier(
    Moment2_Cos_Filter, time2)

Spectrum2_cos_fil = Spectrum2_cos_fil/division_const
Spectrum2_exp_fil = Spectrum2_exp_fil/division_const
Spectrum_c2_cos_fil = Spectrum_c2_cos_fil/division_const
Spectrum_c2_exp_fil = Spectrum_c2_exp_fil/division_const


kind = 'IS'
with open(Path+OutPutFile_IS, 'w') as f:
    ylabel, legend = GetUnits(L_val, kind)
    f.write('{:<10} \t {:<10} \t {:<10}  \n'.format('Energy(MeV)',
            'IS Strength (Exp Filter)', 'IS Strength (Cos Filter)  Units = ['+ylabel[4:]+']'))
    for E, S1, S2 in zip(E1_Exp_fil_IS, Spectrum1_exp_fil_IS, Spectrum1_cos_fil_IS):
        f.write('{:<12.2f} \t {:<18.7f} \t {:<12.7f} \n'.format(E, S1, S2))
kind = 'IV'
with open(Path+OutPutFile_IV, 'w') as f:
    ylabel, legend = GetUnits(L_val, kind)
    f.write('{:<10} \t {:<10} \t {:<10}  \n'.format('Energy(MeV)',
            'IV Strength (Exp Filter)', 'IV Strength (Cos Filter)  Units = ['+ylabel[4:]+']'))
    for E, S1, S2 in zip(E1_Exp_fil_IV, Spectrum1_exp_fil_IV, Spectrum1_cos_fil_IV):
        f.write('{:<12.2f} \t {:<18.7f} \t {:<12.7f} \n'.format(E, S1, S2))

print('\n -------------------------------------------------')
print(
    f'Isoscalar Strength function calculated and written in the file = {Path+OutPutFile_IS}')
print('--------------------------------------------------')
print(
    f'IsoVector Strength function calculated and written in the file = {Path+OutPutFile_IV}')
print('--------------------------------------------------')
