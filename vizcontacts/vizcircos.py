import uuid
import os
import shutil
import math

# write links in main circos conf file
def write_links(conf_file, coupling_dicts_for_sequence_indexed_by_colors, links_folder):

    for color in coupling_dicts_for_sequence_indexed_by_colors:
        links_filename = links_folder+color+".links"
        with open(links_filename, 'w') as f:
            logbase = 5.0
            thick_coeff = 20
            d = coupling_dicts_for_sequence_indexed_by_colors[color]
            for c in d:
                t = 10*d[c]
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
    output_circos_image_tmp = tmp_name+".png"
    write_karyotype(karyotype_filename, sequence)
    write_conf(circos_conf_filename, karyotype_filename, links_folder, output_circos_image_tmp, coupling_dicts_for_sequence_indexed_by_colors)
    os.system("circos -silent -conf "+circos_conf_filename)
    output_circos_image = of+"circos.png"
    shutil.move(output_circos_image_tmp, output_circos_image)
    os.system("xdg-open "+output_circos_image)
