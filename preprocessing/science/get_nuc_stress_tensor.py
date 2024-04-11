### compute nucleation stress tensor 
### created by Jeremy @2024-04

import numpy as np
np.set_printoptions(precision=3)
import argparse
from argparse import RawTextHelpFormatter

def fnomral(dip,strike):
    ### when you look at the strike direction, the fault dips to your right
    dip_rad, strike_rad = np.deg2rad(dip), np.deg2rad(strike)
    n1 = -1*np.sin(dip_rad)*np.sin(strike_rad)
    n2 = np.sin(dip_rad)*np.cos(strike_rad)
    n3 = -1*np.cos(dip_rad)
    return np.array([n1, n2, n3])

def fdip(dip,strike):
    ### when you look at the strike direction, the fault dips to your right
    dip_rad, strike_rad = np.deg2rad(dip), np.deg2rad(strike)
    n1 = -1*np.sin(strike_rad)*np.cos(dip_rad)
    n2 = np.cos(strike_rad)*np.cos(dip_rad)
    n3 = np.sin(dip_rad)
    return np.array([n1, n2, n3])

def fstrike(dip,strike):
    ### when you look at the strike direction, the fault dips to your right
    dip_rad, strike_rad = np.deg2rad(dip), np.deg2rad(strike)
    n1 = np.cos(strike_rad)
    n2 = np.sin(strike_rad)
    n3 = 0
    return np.array([n1, n2, n3])


parser = argparse.ArgumentParser(
    description="Compute nucleation stress tensor given fault geometry and stress level",
    epilog="Examples: \npython get_nuc_stresstensor.py --geometry 9.58 202.27 90  --Traction -18.82 --yaml \npython get_nuc_stresstensor.py --geometry 9.58 202.27 90  --R_nuc 1.1 0.7 0.1 33.6",
    formatter_class=RawTextHelpFormatter
)
parser.add_argument('--geometry',
                    type=float, nargs=3,
                    required=True,
                    help="--geometry dip strike rake", metavar=('DIP', 'STRIKE', 'RAKE'))

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--Traction', type=float, help="Traction [MPa] on the nucleation patch")
group.add_argument('--R_nuc', type=float, nargs=4, metavar=('R', 'MU_S', 'MU_D', 'Pn'),
                   help="--R_nuc R mu_s mu_d Pn (R value, mu_s, mu_d, normal stress [MPa])")

parser.add_argument('--yaml', default=False, action='store_true', help='Generate yaml file for nucleation stress')
args = parser.parse_args()


dip, strike, rake = args.geometry

### compute traction if needed
if args.R_nuc is None:
    T0 = args.Traction
else:
    R_nuc, fs, fd, Pn0 = args.R_nuc
    T0 = R_nuc*(fs - fd)*Pn0 + fd*Pn0

Td0 = T0 * np.sin(np.deg2rad(rake))
Ts0 = T0 * np.cos(np.deg2rad(rake))

### we want S - stress tensor on local coordinates, and backward rotation
S_f = np.array([[0,0,Ts0],[0,0,Td0],[Ts0,Td0,0]])

### fault normal vec and Rotation matrix
nn = fnomral(dip, strike)
nd = fdip(dip, strike)
ns = fstrike(dip, strike)
R = np.vstack([ns,nd,nn])

### we backward rotate S_f from fault local coords to cartesian coords
# foward R@S_cartesian@R.T = S_f 
# backward R.T@S_f@R = S_cartesian
S_cartesian = R.T@S_f@R
print(f'Fault local coordinates (Pn, Td, Ts): 0, {Td0:.2f}, {Ts0:.2f}')
print(f'Projection from Sij: {S_cartesian@nn@nn:.2f}, {S_cartesian@nn@nd:.2f}, {S_cartesian@nn@ns:.2f}, ')

print("Sij =")
print(f" {S_cartesian}")


### generate yaml if needed
yaml_s = f"""!LuaMap
  returns: [nuc_xx, nuc_yy, nuc_zz, nuc_xy, nuc_yz, nuc_xz]
  function: | 
    function f (x)
      xc = -18440.7
      yc = -100744.1
      r =  math.sqrt((x["x"]-xc)^2 +(x["y"]-yc)^2)
      r_crit = 10000.0
      if ( r <= r_crit) then
        ShapeNucleation = math.exp(r^2/(r^2 - r_crit^2))
      else
        ShapeNucleation = 0.0
      end

      return {{
        nuc_xx = ShapeNucleation * {S_cartesian[0,0]:.3f} * 1000000;
        nuc_yy = ShapeNucleation * {S_cartesian[1,1]:.3f} * 1000000;
        nuc_zz = ShapeNucleation * {S_cartesian[2,2]:.3f} * 1000000;
        nuc_xy = ShapeNucleation * {S_cartesian[0,1]:.3f} * 1000000;
        nuc_yz = ShapeNucleation * {S_cartesian[1,2]:.3f} * 1000000;
        nuc_xz = ShapeNucleation * {S_cartesian[0,2]:.3f} * 1000000;
      }}
    end
"""

if args.yaml:
    with open('nucleation_stress.yaml', 'w') as f:
        f.write(yaml_s)
        
    print('wrote nucleation_stress.yaml')



