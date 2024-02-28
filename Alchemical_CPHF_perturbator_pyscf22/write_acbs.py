from AC_basis import p631g,p321g,ccpvDZ
import argparse
import basis_set_exchange as bse


from basis_set_exchange.readers.read  import _reader_map
import numpy as np

valid_formats=list(_reader_map.keys())
valid_formats.remove('genbas')
basis2function={"631g":p631g,
                "321g":p321g,
                "ccpvdz":ccpvDZ
}



def write321g(bs_dict,fz_bs,z):
    #1s
    bs_dict['elements'][z]['electron_shells'][0]['exponents']=[np.format_float_scientific(g[0], precision=12, trim='0') for g in fz_bs[0][1:]]
    bs_dict['elements'][z]['electron_shells'][0]['coefficients']=[[np.format_float_scientific(g[1], precision=12, trim='0') for g in fz_bs[0][1:]]]
    ##2s,2p
    bs_dict['elements'][z]['electron_shells'][1]['exponents']=[np.format_float_scientific(g[0], precision=12, trim='0') for g in fz_bs[1][1:]]
    bs_dict['elements'][z]['electron_shells'][1]['coefficients'][0]=[np.format_float_scientific(g[1], precision=12, trim='0') for g in fz_bs[1][1:]]
    bs_dict['elements'][z]['electron_shells'][1]['coefficients'][1]=[np.format_float_scientific(g[1], precision=12, trim='0') for g in fz_bs[3][1:]]
    ##2s-,2p-
    bs_dict['elements'][z]['electron_shells'][2]['exponents']=[np.format_float_scientific(g[0], precision=12, trim='0') for g in fz_bs[2][1:]]
    bs_dict['elements'][z]['electron_shells'][2]['coefficients'][0]=[np.format_float_scientific(g[1], precision=12, trim='0') for g in fz_bs[2][1:]]
    bs_dict['elements'][z]['electron_shells'][2]['coefficients'][1]=[np.format_float_scientific(g[1], precision=12, trim='0') for g in fz_bs[4][1:]]

    return bs_dict

def write631g(bs_dict,fz_bs,z ):
    #1s
    bs_dict['elements'][z]['electron_shells'][0]['exponents']=[np.format_float_scientific(g[0], precision=10) for g in fz_bs[0][1:]]
    bs_dict['elements'][z]['electron_shells'][0]['coefficients']=[[np.format_float_scientific(g[1], precision=10).replace('e','E') for g in fz_bs[0][1:]]]
    ##2s,2p
    bs_dict['elements'][z]['electron_shells'][1]['exponents']=[np.format_float_scientific(g[0]) for g in fz_bs[1][1:]]
    bs_dict['elements'][z]['electron_shells'][1]['coefficients'][0]=[np.format_float_scientific(g[1]) for g in fz_bs[1][1:]]
    bs_dict['elements'][z]['electron_shells'][1]['coefficients'][1]=[np.format_float_scientific(g[1]) for g in fz_bs[3][1:]]
    ##2s-,2p-
    bs_dict['elements'][z]['electron_shells'][2]['exponents']=[np.format_float_scientific(g[0]) for g in fz_bs[2][1:]]
    bs_dict['elements'][z]['electron_shells'][2]['coefficients'][0]=[np.format_float_scientific(g[1]) for g in fz_bs[2][1:]]
    bs_dict['elements'][z]['electron_shells'][2]['coefficients'][1]=[np.format_float_scientific(g[1]) for g in fz_bs[4][1:]]
    return bs_dict


def write_ccpvdz(bs_dict,fz_bs,z):
    #1s,2s,3s
    bs_dict['elements'][z]['electron_shells'][0]['exponents']=[np.format_float_scientific(g[0], precision=12, trim='0') for g in fz_bs[0][1:]]
    bs_dict['elements'][z]['electron_shells'][0]['coefficients'][0]=[np.format_float_scientific(g[1], precision=12, trim='0') for g in fz_bs[0][1:]]
    bs_dict['elements'][z]['electron_shells'][0]['coefficients'][1]=[np.format_float_scientific(g[2], precision=12, trim='0') for g in fz_bs[0][1:]]
    bs_dict['elements'][z]['electron_shells'][0]['coefficients'][2]=[np.format_float_scientific(g[3], precision=12, trim='0') for g in fz_bs[0][1:]]
    ###2p,3p
    bs_dict['elements'][z]['electron_shells'][1]['exponents']=[np.format_float_scientific(g[0], precision=12, trim='0') for g in fz_bs[1][1:]]
    bs_dict['elements'][z]['electron_shells'][1]['coefficients'][0]=[np.format_float_scientific(g[1], precision=12, trim='0') for g in fz_bs[1][1:]]
    bs_dict['elements'][z]['electron_shells'][1]['coefficients'][1]=[np.format_float_scientific(g[2], precision=12, trim='0') for g in fz_bs[1][1:]]
    ###3d
    bs_dict['elements'][z]['electron_shells'][2]['exponents']=[np.format_float_scientific(g[0], precision=12, trim='0') for g in fz_bs[2][1:]]
    bs_dict['elements'][z]['electron_shells'][2]['coefficients'][0]=[np.format_float_scientific(g[1], precision=12, trim='0') for g in fz_bs[2][1:]]
    return bs_dict

basis2writer={"631g":write631g,
                "321g":write321g,
                "ccpvdz":write_ccpvdz
}



def write_bs(basis,charge,out_format):
    basis_n=basis.lower().replace("-","")
    Z_element=round(charge)
    try:
        basis_f=basis2function[basis_n]
        bs_dict =bse.get_basis(basis,elements=Z_element)
    except:
        print("{} not yet implemented or mistyped. Only available 3-21G, 6-31G,cc-pVDZ".format(basis))
        raise Exception("Invalid basis set")
    fz_bs=basis_f(charge)
    bs_dict=basis2writer[basis_n](bs_dict,fz_bs,str(Z_element))
    extension=_reader_map[out_format]['extension']
    bse.write_formatted_basis_file(bs_dict,f'{basis}_Z_{charge}{extension}',basis_fmt=out_format)
    print(f"The basis file is: {f'{basis}_Z_{charge}{extension}'}")





if __name__=="__main__":
    parser = argparse.ArgumentParser(description="the program writes a basis set file for alchemical consistent basis sets, for a basis set, a nuclear charge and a given output format")

    parser.add_argument("basis", type=str, help="The basis set name")
    parser.add_argument("charge", type=float, help="The nuclear charge")
    parser.add_argument("out_format", type=str, help="The output format")

    args = parser.parse_args()
    if args.charge<3 or args.charge>10 :
        print (f"Basis set for Z= {args.charge } are not yet implemented, 3<Z<10")
        raise Exception("Invalid charge")
    
    if args.out_format not in valid_formats:
        print("BSE supports the following formats: ",valid_formats,"\n")
        print(f"{args.out_format} it is not a valid output format")
        raise Exception("Invalid format")


    write_bs(args.basis, args.charge,args.out_format) 


 

