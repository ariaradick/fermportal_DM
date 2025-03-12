# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 13.3.1 for Mac OS X ARM (64-bit) (July 24, 2023)
# Date: Fri 2 Aug 2024 23:24:09


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot



GC_1 = Coupling(name = 'GC_1',
                value = '-0.3333333333333333*(ee*complex(0,1))',
                order = {'QED':1})

GC_2 = Coupling(name = 'GC_2',
                value = '(-2*ee*complex(0,1))/3.',
                order = {'QED':1})

GC_3 = Coupling(name = 'GC_3',
                value = '(2*ee*complex(0,1))/3.',
                order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = '-(ee*complex(0,1))',
                order = {'QED':1})

GC_5 = Coupling(name = 'GC_5',
                value = 'ee*complex(0,1)',
                order = {'QED':1})

GC_6 = Coupling(name = 'GC_6',
                value = '(8*ee**2*complex(0,1))/9.',
                order = {'QED':2})

GC_7 = Coupling(name = 'GC_7',
                value = 'ee**2*complex(0,1)',
                order = {'QED':2})

GC_8 = Coupling(name = 'GC_8',
                value = '-G',
                order = {'QCD':1})

GC_9 = Coupling(name = 'GC_9',
                value = '-(complex(0,1)*G)',
                order = {'QCD':1})

GC_10 = Coupling(name = 'GC_10',
                 value = 'complex(0,1)*G',
                 order = {'QCD':1})

GC_11 = Coupling(name = 'GC_11',
                 value = '(4*ee*complex(0,1)*G)/3.',
                 order = {'QCD':1,'QED':1})

GC_12 = Coupling(name = 'GC_12',
                 value = 'complex(0,1)*G**2',
                 order = {'QCD':2})

GC_13 = Coupling(name = 'GC_13',
                 value = '-6*complex(0,1)*lam',
                 order = {'QED':2})

GC_14 = Coupling(name = 'GC_14',
                 value = '-(complex(0,1)*lambda1)',
                 order = {'LFDM':1})

GC_15 = Coupling(name = 'GC_15',
                 value = '-(complex(0,1)*lambda2)',
                 order = {'LFDM':1})

GC_16 = Coupling(name = 'GC_16',
                 value = '-(complex(0,1)*lambda3)',
                 order = {'LFDM':1})

GC_17 = Coupling(name = 'GC_17',
                 value = '(ee**2*complex(0,1))/(2.*sw**2)',
                 order = {'QED':2})

GC_18 = Coupling(name = 'GC_18',
                 value = '-((ee**2*complex(0,1))/sw**2)',
                 order = {'QED':2})

GC_19 = Coupling(name = 'GC_19',
                 value = '(cw**2*ee**2*complex(0,1))/sw**2',
                 order = {'QED':2})

GC_20 = Coupling(name = 'GC_20',
                 value = '(ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_21 = Coupling(name = 'GC_21',
                 value = '(cw*ee*complex(0,1))/(2.*sw)',
                 order = {'QED':1})

GC_22 = Coupling(name = 'GC_22',
                 value = '(cw*ee*complex(0,1))/sw',
                 order = {'QED':1})

GC_23 = Coupling(name = 'GC_23',
                 value = '(-2*cw*ee**2*complex(0,1))/sw',
                 order = {'QED':2})

GC_24 = Coupling(name = 'GC_24',
                 value = '-0.16666666666666666*(ee*complex(0,1)*sw)/cw',
                 order = {'QED':1})

GC_25 = Coupling(name = 'GC_25',
                 value = '(ee*complex(0,1)*sw)/(3.*cw)',
                 order = {'QED':1})

GC_26 = Coupling(name = 'GC_26',
                 value = '(2*ee*complex(0,1)*sw)/(3.*cw)',
                 order = {'QED':1})

GC_27 = Coupling(name = 'GC_27',
                 value = '(ee*complex(0,1)*sw)/cw',
                 order = {'QED':1})

GC_28 = Coupling(name = 'GC_28',
                 value = '(-8*ee**2*complex(0,1)*sw)/(9.*cw)',
                 order = {'QED':2})

GC_29 = Coupling(name = 'GC_29',
                 value = '(-4*ee*complex(0,1)*G*sw)/(3.*cw)',
                 order = {'QCD':1,'QED':1})

GC_30 = Coupling(name = 'GC_30',
                 value = '(8*ee**2*complex(0,1)*sw**2)/(9.*cw**2)',
                 order = {'QED':2})

GC_31 = Coupling(name = 'GC_31',
                 value = '-0.5*(cw*ee*complex(0,1))/sw - (ee*complex(0,1)*sw)/(6.*cw)',
                 order = {'QED':1})

GC_32 = Coupling(name = 'GC_32',
                 value = '-0.5*(cw*ee*complex(0,1))/sw + (ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_33 = Coupling(name = 'GC_33',
                 value = '(cw*ee*complex(0,1))/(2.*sw) + (ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_34 = Coupling(name = 'GC_34',
                 value = 'ee**2*complex(0,1) + (cw**2*ee**2*complex(0,1))/(2.*sw**2) + (ee**2*complex(0,1)*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_35 = Coupling(name = 'GC_35',
                 value = '-6*complex(0,1)*lam*vev',
                 order = {'QED':1})

GC_36 = Coupling(name = 'GC_36',
                 value = '(ee**2*complex(0,1)*vev)/(2.*sw**2)',
                 order = {'QED':1})

GC_37 = Coupling(name = 'GC_37',
                 value = 'ee**2*complex(0,1)*vev + (cw**2*ee**2*complex(0,1)*vev)/(2.*sw**2) + (ee**2*complex(0,1)*sw**2*vev)/(2.*cw**2)',
                 order = {'QED':1})

GC_38 = Coupling(name = 'GC_38',
                 value = '-((complex(0,1)*yb)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_39 = Coupling(name = 'GC_39',
                 value = '-((complex(0,1)*yt)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_40 = Coupling(name = 'GC_40',
                 value = '-((complex(0,1)*ytau)/cmath.sqrt(2))',
                 order = {'QED':1})

