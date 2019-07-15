import argparse
import pathlib

from vizpotts.vizpotts import *
from basic_modules.potts_model import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--potts_models', help="Potts model (msgpack files)", type=pathlib.Path, nargs='+')
    #parser.add_argument('-m', '--mode', help="Visualization mode", choices=("whole_potts_model"), default="whole_potts_model")
    parser.add_argument('-1', '--start_at_1', help="Start numbering at 1", action='store_true', default=True), 
    parser.add_argument('-0', '--start_at_0', help="Start numbering at 0", action='store_true', default=False), 
    args = vars(parser.parse_args())


    start_at_1 = args["start_at_1"] and not args["start_at_0"]

    for msgpack in args["potts_models"]:
        mrf = Potts_Model.from_msgpack(msgpack)
        visualize_mrf(mrf, start_at_1=start_at_1, show_figure=False)
    plt.show()



if __name__ == '__main__':
    main()
