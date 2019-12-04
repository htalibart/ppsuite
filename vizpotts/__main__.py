import argparse
import pathlib

from vizpotts.vizpotts import *
from comutils.potts_model import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--potts_models', help="Potts model (msgpack files)", type=pathlib.Path, nargs='+')
    parser.add_argument('-w', '--w_mode', help="Visualize_one_w_matrix", action='store_true', default=False)
    parser.add_argument('-i', '--i_index', help="i index", type=int, default=None)
    parser.add_argument('-j', '--j_index', help="j index", type=int, default=None)
    parser.add_argument('-1', '--start_at_1', help="Start numbering at 1", action='store_true', default=True), 
    parser.add_argument('-0', '--start_at_0', help="Start numbering at 0", action='store_true', default=False), 
    args = vars(parser.parse_args())


    start_at_1 = args["start_at_1"] and not args["start_at_0"]

    if (args["i_index"] is not None) and (args["j_index"] is not None):
        msgpack = args["potts_models"][0]
        mrf = Potts_Model.from_msgpack(msgpack)
        i = args["i_index"]
        j = args["j_index"]
        if start_at_1:
            i-=1
            j-=1
        plot_one_wij(mrf.w[i][j])
    else:
        for msgpack in args["potts_models"]:
            mrf = Potts_Model.from_msgpack(msgpack)
            visualize_mrf(mrf, start_at_1=start_at_1, show_figure=False)
        plt.show()



if __name__ == '__main__':
    main()
