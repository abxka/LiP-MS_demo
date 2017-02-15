load /Users/abdullah/science/uzh/research/public/NatProt/demo/example_output/pdb/4twuA.pdb1
color grey30
hide; set cartoon_fancy_helices, 1; show cartoon; remove solvent; create lig, hetatm; show sticks, lig; show spheres, lig; set sphere_scale, 0.3, lig; cmd.color("grey70", "lig"); util.cnc("lig"); zoom;
hide everything, all and not chain A
bg white
cmd.alias("highlight", "select 1-ALELFR-134-139, resi 134-139 and chain A\ndisable 1-ALELFR-134-139\ncolor yellow, 1-ALELFR-134-139\nlabel 1-ALELFR-134-139 and resi 134 and chain A and name CA, 1\nselect 2-GHHEAELKPL-80-89, resi 80-89 and chain A\ndisable 2-GHHEAELKPL-80-89\ncolor yellow, 2-GHHEAELKPL-80-89\nlabel 2-GHHEAELKPL-80-89 and resi 80 and chain A and name CA, 2\nselect 3-GHHEAELKPLAQ-80-91, resi 80-91 and chain A\ndisable 3-GHHEAELKPLAQ-80-91\ncolor yellow, 3-GHHEAELKPLAQ-80-91\nlabel 3-GHHEAELKPLAQ-80-91 and resi 80 and chain A and name CA, 3\nselect 4-GHHEAELKPLAQSHATK-80-96, resi 80-96 and chain A\ndisable 4-GHHEAELKPLAQSHATK-80-96\ncolor yellow, 4-GHHEAELKPLAQSHATK-80-96\nlabel 4-GHHEAELKPLAQSHATK-80-96 and resi 80 and chain A and name CA, 4\nselect 5-GLSDGEWQQVLNVWGK-1-16, resi 1-16 and chain A\ndisable 5-GLSDGEWQQVLNVWGK-1-16\ncolor yellow, 5-GLSDGEWQQVLNVWGK-1-16\nlabel 5-GLSDGEWQQVLNVWGK-1-16 and resi 1 and chain A and name CA, 5\nselect 6-HGTVVLTALGGILK-64-77, resi 64-77 and chain A\ndisable 6-HGTVVLTALGGILK-64-77\ncolor yellow, 6-HGTVVLTALGGILK-64-77\nlabel 6-HGTVVLTALGGILK-64-77 and resi 64 and chain A and name CA, 6\nselect 7-HPGDFGADAQGAMTK-119-133, resi 119-133 and chain A\ndisable 7-HPGDFGADAQGAMTK-119-133\ncolor yellow, 7-HPGDFGADAQGAMTK-119-133\nlabel 7-HPGDFGADAQGAMTK-119-133 and resi 119 and chain A and name CA, 7\nselect 8-LFTGHPETLEK-32-42, resi 32-42 and chain A\ndisable 8-LFTGHPETLEK-32-42\ncolor yellow, 8-LFTGHPETLEK-32-42\nlabel 8-LFTGHPETLEK-32-42 and resi 32 and chain A and name CA, 8\nselect 9-VEADIAGHGQEVLIR-17-31, resi 17-31 and chain A\ndisable 9-VEADIAGHGQEVLIR-17-31\ncolor yellow, 9-VEADIAGHGQEVLIR-17-31\nlabel 9-VEADIAGHGQEVLIR-17-31 and resi 17 and chain A and name CA, 9\nselect 10-YLEFISDAIIHVLHSK-103-118, resi 103-118 and chain A\ndisable 10-YLEFISDAIIHVLHSK-103-118\ncolor yellow, 10-YLEFISDAIIHVLHSK-103-118\nlabel 10-YLEFISDAIIHVLHSK-103-118 and resi 103 and chain A and name CA, 10\n")
highlight
select nontryptic_is_red, none
select nontryptic_is_red, (resi 89 and chain A) or nontryptic
color red, nontryptic_is_red
select nontryptic_is_red, (resi 91 and chain A) or nontryptic
color red, nontryptic_is_red
select nontryptic_is_red, (resi 1 and chain A) or nontryptic
color red, nontryptic_is_red
select nontryptic_is_red, (resi 89 and chain A) or nontryptic
color red, nontryptic_is_red
select nontryptic_is_red, (resi 91 and chain A) or nontryptic
color red, nontryptic_is_red
select nontryptic_is_red, (resi 1 and chain A) or nontryptic
color red, nontryptic_is_red
select nontryptic_is_red, (resi 89 and chain A) or nontryptic
color red, nontryptic_is_red
select nontryptic_is_red, (resi 91 and chain A) or nontryptic
color red, nontryptic_is_red
select nontryptic_is_red, (resi 1 and chain A) or nontryptic
color red, nontryptic_is_red
disable nontryptic_is_red
