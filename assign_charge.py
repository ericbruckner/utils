def assign_charge(PA_seq, pH = 7, charge_frac = 0.5):
    
    #Extract peptide sequence from PA_sequence
    if (PA_seq[0] == 'C') and PA_seq[1].isdigit():
        if PA_seq[2].isdigit():
            num_alkylC = int(PA_seq[1:3])
            seq = PA_seq[3:]
        else:
            num_alkylC = int(PA_seq[1])
            seq = PA_seq[2:]
    else: # no alkyl
        seq = PA_seq
        num_alkylC = 0
        
    #AA Charge Related Properties
    AA_index = ['A','R','N','D','C','E','Q','G','H','O','I','L','K','M','F','P','U','S','T','W','Y','V']
    AA_pKa = [-1,12.48,-1,3.65,8.18,4.25,-1,-1,6,-1,-1,-1,10.53,-1,-1,-1,-1,-1,-1,-1,10.07,-1] #pKa of Side Chain from D.R. Lide, Handbook of Chemistry and Physics, 72nd Edition, CRC Press, Boca Raton, FL, 1991.
    AA_charge = [0,1,0,-1,0,-1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0]

    #Count Number of Charged Residues
    w = []
    for aa in seq:
        i = AA_index.index(aa)
        if AA_charge[i] == 1:
            if pH > AA_pKa[i]:
                w.append(0)     
            else:
                w.append(AA_charge[i])                          
        elif AA_charge[i] == -1:
            if pH > AA_pKa[i]:
                w.append(AA_charge[i])     
            else:
                w.append(0)
        elif AA_charge[i] == 0:
            w.append(AA_charge[i])       
    charge_count = sum(np.abs(w))

    #Define Charge on AA in self-assembled system
    if pH < 8 and pH > 6:
        if charge_count == 1:
            num_charged_res = np.ceil(charge_count * charge_frac)
        else:
            num_charged_res = np.floor(charge_count * charge_frac)
    else:
        num_charged_res = np.ceil(charge_count * charge_frac)

    num_uncharged_res = charge_count - num_charged_res

    #Build List of Residue Charge
    residue_charge = []
    count = 0
    for i, AA in enumerate(seq):
        if AA in 'KRDE':
            if np.abs(w[i]) == 1 and count < num_uncharged_res:
                residue_charge.append((AA,0))
                count = count + 1
            elif np.abs(w[i]) == 1 and count >= num_uncharged_res:
                residue_charge.append((AA,w[i]))
            else:
                residue_charge.append((AA,0))     

    return residue_charge 
