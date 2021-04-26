# Towards an MD simulation of a crystal of TPA

## Polymorph I (Bailey and Brown, 1967)
Cambridge CCDC: dep no 1269122

Space Group: P 1 (2)

Cell: a 7.730Å b 6.443Å c 3.749Å, α 92.75° β 109.15° γ 95.95° 

## Polymorph II (Bailey and Brown, 1967)
Cambridge CCDC: dep no 1269123

Space Group: P 1 (2)

Cell: a 9.54(1)Å b 5.34(1)Å c 5.02(1)Å, α 86.95(5)° β 134.65(5)° γ 104.90(5)° 

## Polymorph III (Sledz, Janczak, and Kubiak, 2001)

Cambridge CCDC: dep no 154875

Space Group: C 2/m (12)

Cell: a 8.940(2)Å b 10.442(2)Å c 3.7900(10)Å, α 90° β 91.21(3)° γ 90° 

1. Download crystal structure from CCDC -> `154875.cif`
2. OpenBabel: fill unit cell
   - `obabel 154875.cif -opdb -O uc.pdb --fillUC keepconnect`
2. VMD: 
   - `uc.pdb` -> Script `mkucmols.tcl` renames atoms, resid, resnames, removes all H's -> `unitcell_mol1.pdb` and `unitcell_mol2.pdb`
3. OpenBabel to add protons:
   - `obabel unitcell_mol1.pdb -p 1 -O ub7.mol2`
4. CharmmGenFF:
   - paramchem.org: upload `ub7.mol2`, download `ub7.str`
   - manually add H IC's to ub7.str
5. VMD:
   - `ortho_unitcell_mol1.pdb` and `ortho_unitcell_mol2.pdb`-> `mkcrystpsf.tcl` -> `my_ub7_unitcell.psf` and `my_ub7_unitcell.pdb`
   - Test a vacuum MD simulation using `vac.namd` -> `my_ub7_vac.dcd`
   - Script `maksupercellpsf.tcl` makes a supercell -> `my_ub7_supercell.psf`, `my_ub7_supercell.pdb`, and `cell.inp`
   - Test an NPT MD with `solv.namd`
