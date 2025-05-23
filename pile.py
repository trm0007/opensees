import math

import pandas as pd

Depth = [0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
PileDia = [18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18]
FOS = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
SPT = [0, 2, 2, 2, 2, 12, 9, 18, 7, 14, 10, 20, 6, 16, 17, 13, 23, 10, 29, 20, 22]
perimeter_P = [math.pi * N / 12 for N in PileDia]
cross_section_A = [math.pi * N * N * 0.25 / 144 for N in PileDia]
CumDepth = [sum(Depth[:i+1]) for i in range(len(Depth))]
Corrected_SPT = [N if N < 15 else N - (N - 15) / 2 for N in SPT]
AverageSPT = [(Corrected_SPT[i] + Corrected_SPT[i-1]) / 2 if i != 0 and i < len(Corrected_SPT)-0 else Corrected_SPT[i] for i in range(len(Corrected_SPT))]
SurfaceArea_in_Pile = [math.pi * N1 * N2 / 12 for N1, N2 in zip(PileDia, Depth)]
Allowable_Skin_Friction = [0.02 * N1 / N2 for N1, N2 in zip(AverageSPT, FOS)]
Total_Allowable_Skin_Friction = [N1 * N2 for N1, N2 in zip(SurfaceArea_in_Pile, Allowable_Skin_Friction)]
Cum_Skin_Friction = [sum(Total_Allowable_Skin_Friction[:i+1]) for i in range(len(Total_Allowable_Skin_Friction))]
Allowable_End_Bearing_Capacity = [(4 * (N1 * N2 / N3) ) for N1, N2, N3 in zip(Corrected_SPT, cross_section_A, FOS)]
Allowable_working_Load_Ton = [0.75 * (N1 + N2) for N1, N2 in zip(Cum_Skin_Friction, Allowable_End_Bearing_Capacity)]
Allowable_working_Load_Kips = [2.2046226218 * (N1) for N1 in Allowable_working_Load_Ton]
Allowable_bearing_capacity = [N1 / N2 for N1, N2 in zip(Allowable_working_Load_Kips, cross_section_A)]
spacing = 1
SF = 1
Soil_spring_stiffness = [12 * N1 * N2 for N1, N2 in zip(Allowable_bearing_capacity, FOS)]
ks2_3 = [N1 * N2 / 12 for N1, N2 in zip(Soil_spring_stiffness, PileDia)]
ks1 = [N1 * N2 for N1, N2 in zip(Soil_spring_stiffness, perimeter_P)]
end_spring = [N1 * N2 for N1, N2 in zip(Soil_spring_stiffness, cross_section_A)]

data = {
    # "Depth": Depth,
    # "Cum Depth": CumDepth,
    # "S.P.T. Value": SPT,
    # "Corrected S.P.T. Value, N-(N-15)/2": [N if N < 15 else N - (N - 15) / 2 for N in SPT],
    # "Average S.P.T.":  AverageSPT,
    # "Surface Area in Pile (sft)": SurfaceArea_in_Pile,
    # "Allowable Skin Friction, Qa = 0.02N/F.S. (Tsf)": Allowable_Skin_Friction,
    # "Total Allowable Skin Friction, (Ton)": Total_Allowable_Skin_Friction,
    # "Cum Skin Friction, (Ton)": Cum_Skin_Friction,
    # "Allowable End Bearing Capacity, Qa = 4N/F.S.": Allowable_End_Bearing_Capacity,
    # "Allowable working Load (Ton)": Allowable_working_Load_Ton,
    # "Allowable working Load (Kips)": Allowable_working_Load_Kips,
    # "Allowable bearing capacity (Kips/sft)": Allowable_bearing_capacity,
    # "perimeter ,P, ft": perimeter_P,
    # "cross section,A,sft": cross_section_A,
    # "Soil spring stiffness, kip/sft": Soil_spring_stiffness,
    # "ks2,3, kip/ft": ks2_3,
    # "ks1, kip/ft": ks1,
    "end spring, kip": end_spring
}

df = pd.DataFrame(data)
print(df)
df.to_csv("data.csv", index=False)




def momentOfInertia_edge(a, b, d):
    cen_y = 2 * (a * 0.5 * a + b * 0)/(2 * a + b)
    Ixx = 2 * ((d * a **3 /12) + d * a * (a * 0.5 - cen_y) ** 2) + ((b * d **3 /12) + d * b * (cen_y ** 2))
    Iyy = (b * d **3 / 12) + 2 * ((a * d ** 3 / 12) + ((a * d) * (b * 0.5) ** 2))
    return Ixx, Iyy, cen_y

# momentOfInertia_edge(110.35, 140.7, 90.7)

def momentOfInertia_corner(a, b,  d):
    cen_x = (a * 0.5 * a + a * a)/(2 * a + 0)
    cen_y = (a * 0.5 * a + a * a)/(2 * a + 0)
    Ixx = ((d * a **3 /12) + (d * a * (a * 0.5 - cen_x) ** 2)) + ((a * d **3 /12) + d * a * (a * 0.5 - cen_x) ** 2)
    Iyy = ((d * a **3 /12) + (d * a * (a * 0.5 - cen_x) ** 2)) + ((a * d **3 /12) + d * a * (a * 0.5 - cen_x) ** 2)
    return Ixx, Iyy, cen_y

momentOfInertia_corner(110.35, 0,  90.7)

def momentOfInertia_interior(a,b, d):
    cen_x = 0
    cen_y = 0
    Ixx = 2 * (d * a **3 /12)  + 2 * ((a * d **3 /12) + d * a * (a * 0.5 ) ** 2)
    Iyy = 2 * (d * a **3 /12)  + 2 * ((a * d **3 /12) + d * a * (a * 0.5 ) ** 2)
    # print( "cen_y" )
    # print( cen_y )
    # print( "Ixx" )
    # print( Ixx )
    # print( "Iyy" )
    # print( Iyy )
    return Ixx, Iyy, cen_y

def ShearStress_Calculation(footing_width, footing_depth, P, Mxx, Myy,distbetncol_cen_y, dist, Eff_Punching_Perimeter, Avg_Eff_Slab_Thickness, Ixx, Tyy, gama_v1, gama_v2):
    qu = (P * 1000) / (footing_width * footing_depth)
    vu = (P * 1000) - (qu * PunchingArea_length * PunchingArea_width)
    M22 = Mxx - vu * distbetncol_cen_y *0.01
    M33 = Myy
    print("Cm-vu")
    print(vu)
    print( "Cm-Eff_Punching_Perimeter" )
    print( Eff_Punching_Perimeter )
    print( "Cm-Avg_Eff_Slab_Thickness" )
    print( Avg_Eff_Slab_Thickness )
    Vmax = (vu * 1)/(Eff_Punching_Perimeter * Avg_Eff_Slab_Thickness) - ((gama_v1 * M22 * 100 * dist)/Ixx) + ((gama_v2 * M33 * 100 * PunchingArea_length * 0.5)/Iyy)
    return Vmax

# momentOfInertia_interior(140.7,0,  90.7)
# Column Punching Check
Slab_Thickness = 100 #cm
bottom_bar_dia = 18 #mm
Clear_Cover = 7.5 #mm
Col_width = 50 #cm
Col_depth = 50 #cm
footing_width = 700 #cm
footing_depth = 350 #cm
#Calculate distance between the edge of the column and footing
LeftEdge_distance = 50
RightEdge_distance = 50
TopEdge_distance = 85
BottomEdge_distance = 15

fc = 250 #MPa

Avg_Eff_Slab_Thickness = (Slab_Thickness - (Clear_Cover + bottom_bar_dia / 10))  # cm

#Calculate the critical section for punching shear
if LeftEdge_distance < Avg_Eff_Slab_Thickness/2 :
    a1 = LeftEdge_distance
else: a1 = Avg_Eff_Slab_Thickness/2

if RightEdge_distance < Avg_Eff_Slab_Thickness/2 :
    a2 = RightEdge_distance
else: a2 = Avg_Eff_Slab_Thickness/2

if TopEdge_distance < Avg_Eff_Slab_Thickness/2 :
    b1 = TopEdge_distance
else: b1 = Avg_Eff_Slab_Thickness/2

if BottomEdge_distance < Avg_Eff_Slab_Thickness/2 :
    b2 = BottomEdge_distance
else: b2 = Avg_Eff_Slab_Thickness/2
PunchingArea_length = a1 + Col_width + a2
print("cs-PunchingArea_length")
print(PunchingArea_length)

PunchingArea_width = b1 + Col_depth + b2
print("cs-PunchingArea_width")
print(PunchingArea_width)
if a1 < Avg_Eff_Slab_Thickness/2 or a2 < Avg_Eff_Slab_Thickness/2 :
    Eff_Punching_Length = 2 * PunchingArea_length  # cm
else:
    Eff_Punching_Length = PunchingArea_length  # cm
if b1 < Avg_Eff_Slab_Thickness/2 or b2 < Avg_Eff_Slab_Thickness/2:
    Eff_Punching_width = 2 * PunchingArea_width  # cm
else:
    Eff_Punching_width = PunchingArea_width  # cm

Eff_Punching_Perimeter = Eff_Punching_Length + Eff_Punching_width  # cm

Conc_Comp_Strength = 250  # kgf/cm2
Reinforcement_Ratio = 0
Section_Inertia_122 = 57760552.1  # cm*
Section_Inertia_133 = 128870348.2  # cm*
Section_Inertia_123 = 0  # cm*
Gamma = 20.371228
Gamma_v3 = 0.429479
Moment_Mu2 = 72.64  # tonf-m
Moment_Mu3 = -126.9  # tonf-m
Shear_Force = -188.96  # tonf
Unbalanced_Moment_Mu2 = 26.97  # tonf-m
Unbalanced_Moment_Mu3 = -54.5  # tonf-m
Max_Design_Shear_Stress = 12.32  # kgf/cm2
Conc_Shear_Stress_Capacity = 8.32  # kgf/cm2
Punching_Shear_Ratio = 1.48

print("Cs-Eff_Punching_Perimeter")
print(Eff_Punching_Perimeter)


# print("PunchingArea_length")
# print(PunchingArea_length)
# print("PunchingArea_width")
# print(PunchingArea_width)
lmbda = 1
lmbda_s = min(math.sqrt(2 / (1 + (Avg_Eff_Slab_Thickness/2.54) * 0.1)), 1)
beta = max(Col_depth, Col_width) / min(Col_depth, Col_width)
gama_v1 = 1 - (1 / (1 + (2/3)*  math.sqrt(max(PunchingArea_length, PunchingArea_width) / min(PunchingArea_length, PunchingArea_width))))
gama_v2 = 1 - (1 / (1 + (2/3)*  math.sqrt(min(PunchingArea_length, PunchingArea_width) / max(PunchingArea_length, PunchingArea_width))))
print("cs-gama_v1")
print(gama_v1)
print("cs-gama_v2")
print(gama_v2)

if a1 < Avg_Eff_Slab_Thickness/2 or a2 < Avg_Eff_Slab_Thickness/2 :
    if b1 < Avg_Eff_Slab_Thickness / 2 or b2 < Avg_Eff_Slab_Thickness / 2:
        alpha_s = 20
        Ixx, Iyy, cen_y = momentOfInertia_corner(PunchingArea_width, PunchingArea_length, Avg_Eff_Slab_Thickness)
    else:
        alpha_s = 30
        Ixx, Iyy, cen_y = momentOfInertia_corner(PunchingArea_width, PunchingArea_length, Avg_Eff_Slab_Thickness)

elif b1 < Avg_Eff_Slab_Thickness/2 or b2 < Avg_Eff_Slab_Thickness/2 :
    if a1 < Avg_Eff_Slab_Thickness / 2 or a2 < Avg_Eff_Slab_Thickness / 2:
        alpha_s = 20
        Ixx, Iyy, cen_y = momentOfInertia_edge( PunchingArea_width, PunchingArea_length, Avg_Eff_Slab_Thickness )
    else:
        alpha_s = 30
        Ixx, Iyy, cen_y = momentOfInertia_edge(PunchingArea_width, PunchingArea_length, Avg_Eff_Slab_Thickness)

else:
    alpha_s = 40
    Ixx, Iyy, cen_y = momentOfInertia_interior( PunchingArea_width, PunchingArea_length, Avg_Eff_Slab_Thickness )


Vc_min = 0.75 * min(4 * lmbda * lmbda_s * (0.070307 * math.sqrt(250/0.070307)), (2 + 4 / beta) * (4 * lmbda * lmbda_s * (0.070307 * math.sqrt(250/0.070307))), (2 + alpha_s * Avg_Eff_Slab_Thickness / Eff_Punching_Perimeter) * (4 * lmbda * lmbda_s * (0.070307 * math.sqrt(250/0.070307))) )
distbetncol_cen_y = PunchingArea_width - Col_width*0.5 - b2 - cen_y
dist = PunchingArea_width - cen_y
print("Cs-dist")
print(distbetncol_cen_y)
print("Cs-Moment of Inertia, Ixx")
print(Ixx)
print("Cs-Moment of Inertia, Iyy")
print(Iyy)
print("alpha_s")
print(alpha_s)
print("Cs- Vc_min")
print(Vc_min)
P = 208.7 #ton
Mxx = -5.14 #ton-m
Myy = -123.72 #ton-m
Vmax = ShearStress_Calculation(footing_width, footing_depth, P, Mxx, Myy,distbetncol_cen_y, dist, Eff_Punching_Perimeter, Avg_Eff_Slab_Thickness, Ixx, Iyy, gama_v1, gama_v2)
print("Cs-Vmax")
print(Vmax)
