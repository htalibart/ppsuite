import uuid
import os
import shutil
import math

from comutils import files_management as fm
from vizcontacts.contacts_management import *

# write links in main circos conf file
def write_links(conf_file, coupling_dicts_for_sequence_indexed_by_colors, links_folder):

    for color in coupling_dicts_for_sequence_indexed_by_colors:
        links_filename = links_folder+color+".links"
        with open(links_filename, 'w') as f:
            logbase = 5.0
            thick_coeff = 20
            d = coupling_dicts_for_sequence_indexed_by_colors[color]
            for c in d:
                t = 20*d[c]
                if (t>=1):
                    f.write(str(c[0]+1)+" 0 1 "+str(c[1]+1)+" 0 1 thickness="+str(t)+"\n")

            # write link in conf_file
            conf_file.write("""
            <link>
            file          = """+links_filename+"""
            color         = """+color+"""
            radius        = 0.99r
            bezier_radius = 0.1r
            thickness     = 5
            </link>
            """) 



# main circos conf file
def write_conf(circos_conf_filename, karyotype_filename, links_folder, output_circos_image, coupling_dicts_for_sequence_indexed_by_colors):
    with open(circos_conf_filename, 'w') as f:
        f.write("""

        karyotype = """+karyotype_filename+"""

       # <plots>
       # <plot>
       # type = histogram
       # file = histogram_filename
       # r1 = 0.99r
       # r0 = 0.80r
       # thickness = 5
       # color = grey
       # extend_bin = no
       # orientation = in
       # </plot>
       # </plots>

        <ideogram>

        <spacing>
        default = 0.005r
        </spacing>

        # ideogram position, thickness and fill
        radius           = 0.90r
        thickness        = 10p
        fill             = yes


        #stroke_thickness = 1
        #stroke_color     = black

        # ideogram labels
        show_label     = yes
        label_with_tag = yes
        label_font     = light
        label_radius   = dims(ideogram,radius_outer) + 0.05r
        label_center   = yes
        label_size     = 20p
        label_color    = black
        label_parallel = no
        label_case     = upper 


        # ideogram cytogenetic bands, if defined in the karyotype file
        show_bands            = yes
        fill_bands            = yes
        band_stroke_thickness = 2
        band_stroke_color     = grey
        band_transparency     = 1

        </ideogram>

        <image>
        <<include etc/image.conf>>
        file* = """+output_circos_image+"""
        </image>

        <links>""")

        write_links(f, coupling_dicts_for_sequence_indexed_by_colors, links_folder)

        f.write("""
        </links>

        # RGB/HSV color definitions, color lists, location of fonts,
        # fill patterns
        <<include etc/colors_fonts_patterns.conf>> # included from Circos distribution

        # debugging, I/O an dother system parameters
        <<include etc/housekeeping.conf>> # included from Circos distribution

        """
        )


# KARYOTYPE
# chr - id label start end color
def write_karyotype(karyotype_filename, sequence):
    with open(karyotype_filename,"w") as f:
        for pos in range(len(sequence)):
            color = "dgrey"
            f.write("chr - "+str(pos+1)+" "+str(pos+1)+"-"+sequence[pos]+" 0 1 "+color+"\n")


def create_circos(circos_output_folder, coupling_dicts_for_sequence_indexed_by_colors, sequence):
    of = str(circos_output_folder)+'/'
    if not os.path.isdir(of):
        os.mkdir(of)
    karyotype_filename = of+"karyotype.txt"
    circos_conf_filename = of+"main_conf.txt"
    links_folder = of+"links/"
    if not os.path.isdir(links_folder):
        os.mkdir(links_folder)
    tmp_name = str(uuid.uuid4())
    write_karyotype(karyotype_filename, sequence)
    write_conf(circos_conf_filename, karyotype_filename, links_folder, tmp_name+".png", coupling_dicts_for_sequence_indexed_by_colors)
    #os.system("circos -silent -conf "+circos_conf_filename)
    os.system("circos -conf "+circos_conf_filename)
    for extension in [".svg", ".png"]:
        shutil.move(tmp_name+extension, of+"circos"+extension)
    output_circos_image = of+"circos.png"
    os.system("xdg-open "+output_circos_image)


def create_circos_from_comfeature_and_pdb_chain(comfeature, pdb_chain, coupling_sep_min=3, top=20, **args):
    couplings_dict = get_contact_scores_for_sequence(comfeature)
    couplings_dict_with_coupling_sep_min = remove_couplings_too_close(couplings_dict, coupling_sep_min)
    smaller_couplings_dict = OrderedDict({c:couplings_dict_with_coupling_sep_min[c] for c in list(couplings_dict_with_coupling_sep_min)[:top]})
    coupling_dicts_for_sequence_indexed_by_colors = get_colored_true_false_dicts(smaller_couplings_dict, pdb_chain, real_sequence=comfeature.sequence, colors={True:'blue', False:'red'})
    circos_output_folder = str(comfeature.folder.absolute())+"/circos_output"
    create_circos(circos_output_folder, coupling_dicts_for_sequence_indexed_by_colors, comfeature.sequence)



def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--feature_folder', help="Feature folder", type=pathlib.Path)
    parser.add_argument('--pdb_file', help="PDB file", type=pathlib.Path, default=None)
    parser.add_argument('-i', '--pdb_id', help="PDB file")
    parser.add_argument('-cid', '--chain_id', help="PDB chain id", default='A')
    parser.add_argument('-sep', '--coupling_sep_min', help="Min. nb residues between members of a coupling", type=int, default=3)
    parser.add_argument('-n', '--top', help="Nb of couplings displayed", type=int, default=20)
    args = vars(parser.parse_args(args))

    comfeature = ComFeature.from_folder(args['feature_folder'])
    pdb_chain = fm.get_pdb_chain(args['pdb_id'], args['pdb_file'], chain_id=args['chain_id'])
    create_circos_from_comfeature_and_pdb_chain(comfeature, pdb_chain, **args)

if __name__=="__main__":
    main()

