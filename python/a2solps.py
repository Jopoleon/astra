import os, re
import numpy as np
import scipy.interpolate

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def update_line(line, keyword, value, nbc, suffix=""):
    parts = line.split()
    parts[0] = f" {keyword}="

    if nbc == 12:
        if suffix:
            value /= 2
        parts[1] = f"{value}{suffix}, "
        parts[2] = f"{value}{suffix}, "
    else:
        parts[1] = f"{value}{suffix}, "

    return " ".join(parts) + " \n"


def lines2arr(lines_in, aolines, offset, abs_flag=True):
    list_out = []
    for i in range(aolines):
        for word in lines_in[8+i+offset].split():
            val = float(word)
            if abs_flag:
                val = abs(val)
            list_out.append(max(val, 0.01) if not np.isnan(val) else 0.01)
    return list_out


def a2solps(solps_dir, bound, iterate, skips, convect, unique, if_strahl=0):

    f4solps = f'{awd}/a2s.dat'

    dens_bound  = int(bound[0])
    power_bound = int(bound[1])

    print('a2solps')
    with open(f4solps, 'r') as f:
        astraoutput = f.readlines()
    lang  = int(astraoutput[1])
    QE    = np.abs(round(float(astraoutput[2]), 3))
    QI    = np.abs(round(float(astraoutput[3]), 3))
    QN    = np.abs(round(float(astraoutput[4]), 3))
    Qmain = np.abs((float(astraoutput[5])))*1e19
    TE    = np.abs(round(float(astraoutput[6]), 3))
    TI    = np.abs(round(float(astraoutput[7]), 3))

    if (lang % 3) != 0:
        aolines = int(np.ceil(lang/3))
    else:
        aolines = int(lang/3)

    DN  = lines2arr(astraoutput, aolines, 8)
    HE  = lines2arr(astraoutput, aolines, 8+aolines)
    HCI = lines2arr(astraoutput, aolines, 8+2*aolines, abs_flag=False)

    SLAT = float(astraoutput[8 + 3*aolines])
    imps = int(astraoutput[-1])
    posi = 8+3*aolines+1
    df_values = DN
    i = 0
    for i in range(imps):
        for j in range(aolines):
            for word in astraoutput[posi + i*aolines + j].split():
                df_num = float(word)
                if df_num == 0:
                    df_num = 0.01
                df_values.append(df_num)
    df = np.reshape(df_values, (lang, imps+1), order='F')
    if imps == 0:
        new_posi = posi + i*aolines
    else:
        new_posi = posi + i*aolines + j + 1
    qf_values = [Qmain]

    for i in range(imps):
        qf_values.append(((float(astraoutput[new_posi+i])))*1e19)
    if imps == 0 :
        rho_pos = new_posi + i
    else:
        rho_pos = new_posi + i + 1
    astra_rhos = []
    for j in range(aolines):
        for word in astraoutput[rho_pos+j].split()
            astra_rhos.append(float(word))

    dens_pos = rho_pos + j + 1
    Nmain = float(astraoutput[dens_pos])*1e19
    N_values = [Nmain]
    for i in range(imps):
        N_values.append(float(astraoutput[dens_pos+1+i])*1e19)
    z_pos = dens_pos + 1 + i
    z_values = [1]
    for i in range(imps):
        z_values.append(int(astraoutput[z_pos+i+1]))
    con_pos = z_pos + 2 + i
    con_values = []
    for i in range(imps+1):
        for j in range(aolines):
            for word in astraoutput[con_pos+i*aolines+j].split()
                con_values.append(float(word))
    con = np.reshape(con_values, (lang, imps+1), order='F')

    XUPAR = HCI

    os.chdir(solps_dir)

    with open('b2mn.dat', 'r') as f:
        mnlines = f.readlines()

    tausolps = None
    has_b2tqna = False

    for i, line in enumerate(mnlines):
        if 'b2mndr_dtim' in line:
            tausolps = float(line.split()[1].strip("'"))
        if 'b2tqna_inputfile' in line:
            has_b2tqna = True
            if '0' in line:
                mnlines[i] = line.replace('0', '1')

    tauastra = 1.0e-3
    nsolps = int(iterate)
    if nsolps == 0:
        nsolps = int((skips * tauastra) / tausolps)

    if not has_b2tqna:
        mnlines.append("'b2tqna_inputfile'          '1'  # Look for b2.transport.inputfile\n")
    mnlines.append(f"b2mndr_ntim'                      '{nsolps}'         # nsteps, if 0 then only write state file \n")

    with open('b2mn.dat', 'w') as f:
        f.writelines(mnlines)

    with open("b2.boundary.parameters") as f:
        boundary_file = f.readlines()

    nbc = int(boundary_file[1].replace(",", "").split()[-1])

# Select values according to power_bound
    if power_bound == 0:
        electron_value = QE
        ion_value = QI
        suffix = "E+06"
    else:
        electron_value = TE
        ion_value = TI
        suffix = ""

    for i, line in enumerate(boundary_file):
        key = line.lower().replace(" ", "")

        if "enepar(1,1)" in key:
            boundary_file[i] = update_line(line, "enepar(1, 1)", electron_value, nbc, suffix)

        elif "enipar(1,1)" in key:
            boundary_file[i] = update_line(line, "enipar(1, 1)", ion_value, nbc, suffix)

    z_count = 0
    # Sparses through the b2fstati and gets the charges
    fb2f = open('b2fstati', 'r')
    for line in fb2f:
        if 'zamin' in line:
            zas = []
            while z_count < 1:
                lines = (fb2f.readline())
                if 'zamax' in lines:
                    break
                else:
                    zas.append(lines)
    fb2f.close()

    zamins = [int(float(value)) for line in zas for value in line.split()]
    neuts = [i for i, z in enumerate(zamins) if z == 0]

    zerostr = ' 0.00, '
    conparlist  = [' conpar(0, 1, 1) =']
    conpar2list = [' conpar(0, 2, 1) =']
    bcconlist   = [' bccon(0, 1) =']
    bccon2list  = [' bccon(0, 2) =']
    conzero = ' 0, '
    con2    = ' 2, '

    for zamin in zamins:
        conparlist.append(zerostr)
        conpar2list.append(zerostr)
        if zamin == 0:
            bcconlist.append(conzero)
            bccon2list.append(conzero)
        else:
            bcconlist.append(con2)
            bccon2list.append(con2)
    for i, neut in enumerate(neuts):
        idx = neut + 1 + z_values[i]
        if dens_bound == 0:
            if i > 0 and unique == 1:
                value = N_values[i]
                bc = 1
            else:
                value = qf_values[i]
                bc = 8
        else:
            if i > 0 and unique == 1:
                value = qf_values[i]
                bc = 8
            else:
                value = N_values[i]
                bc = 1
        conparlist[idx]  = f' {value:.3E}, '
        conpar2list[idx] = f' {value:.3E}, '
        bcconlist[idx]   = f' {bc}, '
        bccon2list[idx]  = f' {bc}, '
    conparlist [-1] += ' \n'
    conpar2list[-1] += ' \n'
    bcconlist  [-1] += ' \n'
    bccon2list [-1] += ' \n'
        
    conparstr  = ' '.join(conparlist)
    conpar2str = ' '.join(conpar2list)
    bcconstr   = ' '.join(bcconlist)
    bccon2str  = ' '.join(bccon2list)
    for x in range(len(boundary_file)):
        if 'conpar(0, 1, 1)' in boundary_file[x].lower().replace(" ", ""):
            boundary_file[x] = conparstr
        if 'bccon(0, 1)' in boundary_file[x].lower().replace(" ", ""):
            boundary_file[x] = bcconstr
        if nbc == 12:
            if 'conpar(0, 2, 1)' in boundary_file[x].lower().replace(" ", ""):
                boundary_file[x] = conpar2str
            if 'bccon(0, 2)' in boundary_file[x].lower().replace(" ", ""):
                boundary_file[x] = bccon2str

    with open('b2.boundary.parameters', 'w') as f:
        f.writelines(boundary_file)

    rho_nums = int(f.readline())
    solps_rhos = []
    with open('SOLPS_rho.txt', 'r') as f:
        for line in f:
            split = line.split()
            for i in range(len(split)):
                solps_rhos.append(float(split[i]))

    cuts = {}

    with open("b2fgmtry") as f:
        for line in f:
            if "nncut" in line:
                nncut = int(next(f))
            elif "leftcut" in line:
                n = int(line.split()[-2])
                cuts["left"] = list(map(int, next(f).split()))[:n]
            elif "rightcut" in line:
                n = int(line.split()[-2])
                cuts["right"] = list(map(int, next(f).split()))[:n]
            elif "topcut" in line:
                n = int(line.split()[-2])
                cuts["top"] = list(map(int, next(f).split()))[:n]

    if nncut == 1:
        core_left  = cuts["left" ][0] + 1
        core_right = cuts["right"][0] + 1
        core_top   = cuts["top"  ][0] + 1
    elif nncut == 2:
        core_left_left   = cuts["left"][0] + 1
        core_left_right  = cuts["left"][1] + 1
        core_right_left  = min(cuts["right"]) + 1
        core_right_right = max(cuts["right"]) + 1
        core_top         = min(cuts["top"]) + 1

    solps_pos = np.loadtxt('SOLPS_r.txt')
    astra_r = scipy.interpolate.interp1d(solps_rhos, solps_pos[:, 0], 'linear', fill_value='extrapolate')(astra_rhos)

    rsep = astra_r[-1]

    r_rsep = astra_r - rsep
    r_rsep = -np.abs(r_rsep)

    solps_r_rsep = scipy.interpolate.interp1d(astra_r, r_rsep, 'linear', fill_value='extrapolate')(solps_pos[:core_top, 0])
    solps_r_rsep = np.round(solps_r_rsep, 5)
    og_solps_r_rsep = solps_r_rsep[0]

    HCI   = scipy.interpolate.interp1d(r_rsep, HCI  , 'linear', fill_value='extrapolate')(solps_r_rsep)
    HE    = scipy.interpolate.interp1d(r_rsep, HE   , 'linear', fill_value='extrapolate')(solps_r_rsep)
    XUPAR = scipy.interpolate.interp1d(r_rsep, XUPAR, 'linear', fill_value='extrapolate')(solps_r_rsep)

    ndata_index = []
    ndata_lines = []

    # Key: (3rd_index, 4th_index) -> list of (positional_value, final_value)
    # Only stores entries where positional_value > 0
    tdata_sol_map = {}

    try:
        with open('b2.transport.inputfile', 'r') as f:
            transport_lines = f.readlines()

            for line in transport_lines:
                if 'tdata' not in line:
                    if 'ndata' in line:
                        ndata_index.append(transport_lines.index(line))
                        ndata_lines.append(line)
                    continue

                # Each line has two tdata entries:
                # tdata(1, seq, idx3, idx4) = positional_value , tdata(2, seq, idx3, idx4) = final_value
                # Find all tdata entries: grab index tuple and value
                # Pattern: tdata(a, b, c, d) = value
                entries = re.findall(
                    r'tdata\s*\(\s*(\d+)\s*, \s*(\d+)\s*, \s*(\d+)\s*, \s*(\d+)\s*\)\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', 
                    line
                )

                if not entries:
                    continue

                # entries is a list of (a, b, c, d, value) tuples
                # First entry (a=1) holds the positional value
                # Second entry (a=2) holds the final/data value
                # Group them by sequence index b — both entries on the line share same b, c, d
                entry_dict = {}
                for (a, b, c, d, val) in entries:
                    a, b, c, d = int(a), int(b), int(c), int(d)
                    val = float(val)
                    key_seq = (b, c, d)
                    if key_seq not in entry_dict:
                        entry_dict[key_seq] = {}
                    entry_dict[key_seq][a] = val

                for (b, c, d), vals in entry_dict.items():
                    if 1 in vals and 2 in vals:
                        posi_val  = vals[1]   # positional value (a=1)
                        final_val = vals[2]   # data value (a=2)
                        key = (c, d)          # (3rd index, 4th index)

                        if posi_val > 0.0001:
                            if key not in tdata_sol_map:
                                tdata_sol_map[key] = []
                            tdata_sol_map[key].append((posi_val, final_val))

    except IOError as e:
        print("Error:", e)

    # ndata_sol_size per (3rd, 4th) combination
    ndata_sol_size = {key: len(entries) for key, entries in tdata_sol_map.items()}
  
    print("ndata_sol_size:", ndata_sol_size)
    def get_sol_size(idx3, idx4):
        return ndata_sol_size.get((idx3, idx4), 0)

    def get_sol_entries(idx3, idx4):
        return tdata_sol_map.get((idx3, idx4), [])

    def write_tdata_block(transport_file, idx2, idx4, neg_positions, neg_values):
        """
        Writes an ndata + tdata block.
        idx2          : second/third index in ndata/tdata (same value, e.g. 1, 3, 4, 6, 7)
        idx4          : fourth index (k)
        neg_positions : list of positional values <= 0 (your new data)
        neg_values    : list of corresponding data values
        """
        sol_entries   = get_sol_entries(idx2, idx4)
        sol_size      = get_sol_size(idx2, idx4)
        neg_size      = len(neg_positions)
        combined_size = neg_size + sol_size

        transport_file.append(f" ndata(1, {idx2} , {idx4} )= {combined_size} , \n")

        # Write new negative/zero position data
        for i in range(neg_size):
            transport_file.append(
                f" tdata(1, {i+1} , {idx2} , {idx4} )= {neg_positions[i]} , "
                f" tdata(2, {i+1} , {idx2} , {idx4} )= {neg_values[i]} , \n"
            )

        # Write positive position data preserved from file
        for w, (posi_val, final_val) in enumerate(sol_entries, start=neg_size + 1):
            transport_file.append(
                f" tdata(1, {w} , {idx2} , {idx4} )= {posi_val} , "
                f" tdata(2, {w} , {idx2} , {idx4} )= {final_val} , \n"
            )

    transport_file = list()
    transport_file.append('&TRANSPORT\n')
    print(solps_r_rsep)
    # --- Index 1 (diffusion) ---
    df_counter = -1
    for k in range(len(zamins)):
        if zamins[k] == 0:
            df_counter += 1
            solps_r_rsep[0] = og_solps_r_rsep
            df0 = scipy.interpolate.interp1d(r_rsep, df[:, df_counter], 'linear', fill_value='extrapolate')(solps_r_rsep)
            solps_r_rsep[0] = -0.4
        else:
            neg_mask = [r <= 0.0001 for r in solps_r_rsep]
            neg_pos  = [solps_r_rsep[i] for i in range(len(solps_r_rsep)) if neg_mask[i]]
            neg_vals = [df0[i] for i in range(len(solps_r_rsep)) if neg_mask[i]]
            write_tdata_block(transport_file, 1, k, neg_pos, neg_vals)

    # --- Index 3 (HCI) ---
    for k in range(len(zamins)):
        if zamins[k] != 0:
            neg_mask = [r <= 0.0001 for r in solps_r_rsep]
            neg_pos  = [solps_r_rsep[i] for i in range(len(solps_r_rsep)) if neg_mask[i]]
            neg_vals = [HCI[i] for i in range(len(solps_r_rsep)) if neg_mask[i]]
            write_tdata_block(transport_file, 3, k, neg_pos, neg_vals)

    # --- Index 4 (HE) ---
    neg_mask = [r <= 0.0001 for r in solps_r_rsep]
    neg_pos  = [solps_r_rsep[i] for i in range(len(solps_r_rsep)) if neg_mask[i]]
    neg_vals = [HE[i] for i in range(len(solps_r_rsep)) if neg_mask[i]]
    write_tdata_block(transport_file, 4, 1, neg_pos, neg_vals)

    # --- Index 6 (convective) ---
    if convect == 1:
        con_counter = -1
        for k in range(len(zamins)):
            if zamins[k] == 0:
                con_counter += 1
                solps_r_rsep[0] = og_solps_r_rsep
                con0 = scipy.interpolate.interp1d(r_rsep, con[:, con_counter], 'linear', fill_value='extrapolate')(solps_r_rsep)
                solps_r_rsep[0] = -0.4
            else:
                neg_mask = [r <= 0.0001 for r in solps_r_rsep]
                neg_pos  = [solps_r_rsep[i] for i in range(len(solps_r_rsep)) if neg_mask[i]]
                neg_vals = [con0[i] for i in range(len(solps_r_rsep)) if neg_mask[i]]

                if get_sol_size(6, k) > 0:
                    write_tdata_block(transport_file, 6, k, neg_pos, neg_vals)
                else:
                    # No existing data in file — write with just ASTRA data
                    transport_file.append(f" ndata(1, 6 , {k} )= {len(neg_pos)} , \n")
                    for i in range(len(neg_pos)):
                        transport_file.append(
                            f" tdata(1, {i+1} , 6 , {k} )= {neg_pos[i]} , "
                            f" tdata(2, {i+1} , 6 , {k} )= {neg_vals[i]} , \n"
                        )

    # --- Index 7 (XUPAR) ---
    if convect == 1:
        for k in range(len(zamins)):
            if zamins[k] != 0:
                neg_mask = [r <= 0.0001 for r in solps_r_rsep]
                neg_pos  = [solps_r_rsep[i] for i in range(len(solps_r_rsep)) if neg_mask[i]]
                neg_vals = [XUPAR[i] for i in range(len(solps_r_rsep)) if neg_mask[i]]
                if get_sol_size(7, k) > 0:
                    write_tdata_block(transport_file, 7, k, neg_pos, neg_vals)
                else:
                    # No existing data in file — write with just ASTRA data
                    transport_file.append(f" ndata(1, 7 , {k} )= {len(neg_pos)} , \n")
                    for i in range(len(neg_pos)):
                        transport_file.append(
                            f" tdata(1, {i+1} , 7 , {k} )= {neg_pos[i]} , "
                            f" tdata(2, {i+1} , 7 , {k} )= {neg_vals[i]} , \n"
                        )

    transport_file.append(' no_pflux =.true. , \n')
    transport_file.append(' /\n')

    # --- Write to file ---
    with open('b2.transport.inputfile', 'w') as f:
        f.writelines(transport_file)
    return
