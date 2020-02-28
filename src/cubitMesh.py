import os
import sys
import argparse

def main(args):
    pathToTrelis = r"C:\Program Files\Cubit 15.5\bin"
    sys.path.append(pathToTrelis)
    # Initialize Cubit
    import cubit
    cubit.init(['cubit', '-nojournal', '-noecho','-nobanner','-nographics', '-warning', 'off', '-information','off'])
    cubit.cmd('reset')
    # Unpack arguments
    xmin, xmax, ymin, ymax = [args.xmin, args.xmax, args.ymin, args.ymax]
    NX, NY = [args.nx, args.ny]
    degree = args.degree
    outname = args.outname
    width = xmax - xmin
    height = ymax - ymin
    # Create surface
    cubit.cmd('create surface rectangle width ' + str(width) + ' height ' + str(height) +' zplane ')
    cubit.cmd("rotate Surface 1  angle 180  about Z include_merged")
    cubit.cmd('move vertex ( at ' + str(-width/2) + ' ' + str(-height/2) + ' 0 ordinal 1 ) location ' + str(xmin) + ' ' + str(ymin) + ' 0 include_merged ')
    # Create block
    cubit.cmd('block 1 surf 1')
    cubit.cmd('block 1 name ' + '"body"')
    cubit.cmd('block 1 element type QUAD4')
    # Create sideset
    cubit.cmd('sideset 1 add curve 2')
    cubit.cmd('sideset 2 add curve 4')
    cubit.cmd('sideset 3 add curve 3')
    cubit.cmd('sideset 4 add curve 1')
    cubit.cmd('sideset 1 name ' + '"Left"')
    cubit.cmd('sideset 2 name ' + '"Right"')
    cubit.cmd('sideset 3 name ' + '"Bottom"')
    cubit.cmd('sideset 4 name ' + '"Top"')
    # Create nodeset
    cubit.cmd('nodeset 1 add curve 2  ')
    cubit.cmd('nodeset 2 add curve 4  ')
    cubit.cmd('nodeset 3 add curve 3  ')
    cubit.cmd('nodeset 4 add curve 1  ')
    cubit.cmd('nodeset 1 name ' + '"Left"')
    cubit.cmd('nodeset 2 name ' + '"Right"')
    cubit.cmd('nodeset 3 name ' + '"Bottom"')
    cubit.cmd('nodeset 4 name ' + '"Top"')
    # Specify mesh density - NX
    cubit.cmd('curve all in sideset 3 4 interval ' + str(NX))
    cubit.cmd('curve all in sideset 3 4  scheme equal')
    cubit.cmd('mesh curve all in sideset 3 4')
    # Specify mesh density - NY
    cubit.cmd('curve all in sideset 1 2 interval ' + str(NY))
    cubit.cmd('curve all in sideset 1 2  scheme equal')
    cubit.cmd('mesh curve all in sideset 1 2')
    # Mesh surface
    cubit.cmd('mesh surf 1')
    # Export Genesis
    cubit.cmd('set exodus netcdf4 on')
    cubit.cmd('export mesh "./' + outname + '.g"  dimension 2  overwrite ')

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xmin",    default=0.,     type=float)
    parser.add_argument("--xmax",    default=1.,     type=float)
    parser.add_argument("--ymin",    default=0.,     type=float)
    parser.add_argument("--ymax",    default=1.,     type=float)
    parser.add_argument("--nx",      default=10,     type=int)
    parser.add_argument("--ny",      default=10,     type=int)
    parser.add_argument("--degree",  default=1,      type=int)
    parser.add_argument("--outname", default="mesh", type=str)
    return parser

if __name__ == "__main__":
    parser = parseArguments()
    args = parser.parse_args()
    main(args)