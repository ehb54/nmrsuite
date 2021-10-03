from os.path import join
raw_dir = "/opt/genapp/rotdif/output/html5/results/users/ychen151/no_project_specified/aug4_ELM"
with open(join(raw_dir,"str_abc.txt")) as in_f:
    sa = in_f.readline()
    fsa = [float(s) for s in sa.split(' ')]
print(type(fsa[0]))
cmx, cmy, cmz, a1, a2, a3, deg1, deg2, deg3 = 1,2,3,4,5,6,7,8,9
to_append = """
#1. Center of Mass of the molecule
cmx, cmy, cmz = %0.3f, %0.3f, %0.3f
#2. Ellipsoid semiaxes length (in Angstrom)
a1,a2,a3 = %0.3f, %0.3f, %0.3f
#3. Color: see https://pymolwiki.org/index.php/Color_Values for more options
color = [0.85, 0.85, 1.00]
#4. Rotation Input: three Euler angles: alpha, beta, gamma (in degrees)
rotationInput = [%0.3f, %0.3f, %0.3f]
tmp = drawEllipsoid(color, cmx, cmy, cmz, a1, a2, a3, *rotationMatrix(rotationInput))
cmd.load_cgo(tmp, 'ellipsoid-cgo')
cmd.set('cgo_transparency', 0.5, 'ellipsoid-cgo')"""%(cmx, cmy, cmz, a1, a2, a3, deg1, deg2, deg3)
with open('test_encode.txt','w') as out_f:
    out_f.write(to_append)
