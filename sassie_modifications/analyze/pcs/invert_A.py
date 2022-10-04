import numpy as np

def invert_A(A=None):

    con_num = np.linalg.cond(A,p=None)# Condition number

    [u, w, v] = np.linalg.svd(A, full_matrices=True, compute_uv=True)# z=u.w.transpose(v)

    print(np.shape(u))
    new_w = np.full((55, 5), 0.0)
    for i in range(2): new_w[i][i] = w[i]
    w = new_w
    # numpy returns a vector, matlab returns a diagonalized matrix
    v = np.transpose(v) # numpy returns the transpose of what matlab returns
    print(u)
    print(w)
    print(v)

    max_w = np.max(w)

    s = np.full((5, 5), 0.0)

    tol = 0.01 * max_w

    for i in range(2): #for loop not inclusive of 5
        if w[i][i] > tol:
            s[i][i] = 1 / (w[i][i])    # Inverse of singular values
        else:
            s[i][i]= 0.0        # Exclude very small values

    A_inv = np.matmul(v, s)# v.s.transpose(u)
    A_inv = np.dot(A_inv, np.transpose(u))

    '''
    print(s)
    print(np.transpose(u))
    print (s * np.transpose(u,axes=None))
    print(np.matmul(s, np.transpose(u)))
    print (A_inv)
    print (w)
    print (v)
    print (con_num)
    '''

    return [A_inv,w,v,con_num]

    #=======================================================================
'''
a = np.array([[0.9614103 , 0.90455435, 0.99458217, 0.31737855, 0.45671415],
       [0.12040298, 0.84392798, 0.76281861, 0.25349974, 0.97096536],
       [0.97279915, 0.27278921, 0.43858012, 0.05284235, 0.09809111],
       [0.33026723, 0.16024311, 0.29228199, 0.56830532, 0.23192066],
       [0.35259912, 0.69101618, 0.18089714, 0.77224122, 0.69253435],
       [0.77895636, 0.78647306, 0.76349102, 0.70765827, 0.50017265],
       [0.40588817, 0.85204723, 0.91226632, 0.3842231 , 0.21995783],
       [0.11646815, 0.3064179 , 0.69549178, 0.15915652, 0.05393551],
       [0.96318775, 0.46032371, 0.83506767, 0.53555313, 0.01730094],
       [0.61638335, 0.02688136, 0.0369607 , 0.97012081, 0.18417634],
       [0.02175112, 0.16465837, 0.63577357, 0.35041669, 0.25288747],
       [0.1488957 , 0.74792589, 0.86378564, 0.55920303, 0.96119228],
       [0.11060531, 0.57356521, 0.5450523 , 0.75630361, 0.4224476 ],
       [0.51202017, 0.35943428, 0.22651565, 0.20780728, 0.58559785],
       [0.81979152, 0.64605513, 0.7469593 , 0.07894681, 0.18905192],
       [0.57384634, 0.0063519 , 0.92023696, 0.57142295, 0.5695416 ],
       [0.94062953, 0.62151824, 0.93731288, 0.41448092, 0.13940826],
       [0.74194169, 0.49688075, 0.13095444, 0.21770585, 0.82176905],
       [0.24400168, 0.44020594, 0.63290417, 0.59041577, 0.09778248],
       [0.8604161 , 0.09714069, 0.37697765, 0.69142945, 0.21966259],
       [0.46220002, 0.34921751, 0.95736011, 0.02068288, 0.78625365],
       [0.10558957, 0.23092473, 0.83417352, 0.0813662 , 0.32764497],
       [0.42142224, 0.56015884, 0.20158301, 0.35181316, 0.47993089],
       [0.84913992, 0.90882509, 0.80503717, 0.27872725, 0.94979299],
       [0.41723691, 0.45107572, 0.95270153, 0.77786851, 0.66862611],
       [0.34880992, 0.857237  , 0.63695915, 0.77469074, 0.14781722],
       [0.39277179, 0.35944608, 0.87984919, 0.14171074, 0.49163577],
       [0.2500409 , 0.44790555, 0.30882482, 0.25815118, 0.23817293],
       [0.84678098, 0.83664   , 0.22335713, 0.01247071, 0.51396654],
       [0.15427164, 0.41216833, 0.46111806, 0.1874917 , 0.38424292],
       [0.89567377, 0.55030944, 0.61160675, 0.37753013, 0.26837536],
       [0.36921814, 0.88229178, 0.77243492, 0.04192881, 0.94692699],
       [0.75891995, 0.73126445, 0.45372556, 0.74074079, 0.82161232],
       [0.03932268, 0.39576279, 0.80276919, 0.23715146, 0.10198288],
       [0.97213167, 0.80906756, 0.34132368, 0.48023556, 0.60360861],
       [0.51833069, 0.93179444, 0.5989958 , 0.36765403, 0.55307576],
       [0.83188919, 0.14480363, 0.46297717, 0.68816116, 0.61226818],
       [0.08188345, 0.79242268, 0.28728244, 0.62496106, 0.10306321],
       [0.30511795, 0.35236937, 0.22695933, 0.69545178, 0.32293616],
       [0.92055052, 0.77417435, 0.02586332, 0.18240288, 0.87008196],
       [0.67107526, 0.1159119 , 0.05807383, 0.77148419, 0.36362567],
       [0.48557226, 0.75524861, 0.71540167, 0.10193888, 0.56657575],
       [0.83855374, 0.06941731, 0.06955075, 0.51837168, 0.21495469],
       [0.84050397, 0.91141171, 0.09858582, 0.36825397, 0.20888214],
       [0.02466058, 0.01636702, 0.65533639, 0.51478374, 0.29634713],
       [0.25049117, 0.27205837, 0.64927188, 0.7250207 , 0.89233049],
       [0.84170616, 0.33509093, 0.68446772, 0.95656012, 0.12452294],
       [0.67474598, 0.15463275, 0.42290742, 0.59249471, 0.03077629],
       [0.30549386, 0.22574406, 0.26957699, 0.16830233, 0.52316017],
       [0.58260703, 0.76268924, 0.20291167, 0.31481827, 0.72707428],
       [0.67036109, 0.27427542, 0.08089001, 0.74531987, 0.40551471],
       [0.82307068, 0.78534842, 0.25073847, 0.42958211, 0.91387073],
       [0.18395594, 0.70507331, 0.93854945, 0.30591785, 0.30307243],
       [0.72918244, 0.98012839, 0.0443729 , 0.27634786, 0.95135521],
       [0.77643436, 0.15877296, 0.50619677, 0.80016206, 0.06064017]])
b = np.array([[1, 2], [3, 4]])
invert_A(b)
'''
