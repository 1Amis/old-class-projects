(* Keplerian orbits *) 

(*plotting a radial orbit*)

Clear["Global`*"];  x0 = 10 ; v0 = 1/10 ; 
h = x0*v0 ; p = h^2 ;  e = Sqrt[1 - p/a] ; 
E0 = (1/2)*(v0)^2 - 1/x0 ;  a = -1/(2 E0);
P = 2*Pi*a^(3/2) ;  n = 300;
(* defne r\[LongEqual]p/(1-\[Epsilon]cos\[Theta]),  and to plot in \
cartesian, use X & Y parametric plot with   
 X\[LongEqual](r)*cos\[Theta],Y\[LongEqual]r*sin\[Theta] ;  *)
xp[i_] = p*Cos[i*2*Pi/n]/(1 - e*Cos[i*2*Pi/n]) ; 
yp[i_] = p*Sin[i*2*Pi/n]/(1 - e*Cos[i*2*Pi/n]) ; 
fullplot = Table[{xp[i], yp[i]}, {i, 1, 300}]; 
ListPlot[fullplot, PlotMarkers -> " 〇 "]

(*numerically solving a diffeq to study the Kepler problem in 2D*)

(*Cromer Algorithm*)

x0 = 10 ; v0 = 0.1 ; 
E0 = N[(1/2)*(v0)^2 - 1/x0 ]; a = N[-1/(2*E0)];
P = N[2*Pi*a^(3/2)] ; dt = P/500; 
rx[1] = 10; ry[1] = 0; vx[1] = 0; vy[1] = 0.1; 
Csolutionr = {{rx[1], ry[1]}};
Csolutionv = {{vx[1], vy[1]}}; n = 1;
Do[
  nn = n + 1;
  {vx[nn], 
    vy[nn]} = {vx[n], vy[n]} - {rx[n], ry[n]}*
     dt/(rx[n]^2 + ry[n]^2)^(3/2) ; 
  {rx[nn], ry[nn]} = {rx[n], ry[n]} + {vx[nn], vy[nn]}*dt; 
  Csolutionr = Append[Csolutionr, {rx[nn], ry[nn]}];
  Csolutionv = Append[Csolutionv, {vx[nn], vy[nn]}];
  {vx[n], vy[n]} = {vx[nn], vy[nn]};
  {rx[n], ry[n]} = {rx[nn], ry[nn]}; , {n, 2500}];
Dimensions[Csolutionr];
ListPlot[{fullplot , Csolutionr},  {PlotStyle -> {Blue, Black}, 
  PlotMarkers -> {〇, None},
  Joined -> {False, True}, PlotLabel -> "Cromer Method", 
  ImageSize -> Large} ]

(*Newton's Method*)

x0 = 10 ; v0 = 0.1 ; 
E0 = N[(1/2)*(v0)^2 - 1/x0 ]; a = N[-1/(2*E0)];
P = N[2*Pi*a^(3/2)] ; dt = N[P/500]; 
rx[1] = 10; ry[1] = 0; vx[1] = 0; vy[1] = 0.1; 
Nsolutionr = {{rx[1], ry[1]}};
Nsolutionv = {{vx[1], vy[1]}}; n = 1;
Do[
  nn = n + 1;
  {rx[nn], ry[nn]} = {rx[n], ry[n]} + {vx[n], vy[n]}*dt; 
  {vx[nn], 
    vy[nn]} = {vx[n], vy[n]} - {rx[nn], ry[nn]}*
     dt/(rx[nn]^2 + ry[nn]^2)^(3/2) ; 
  Nsolutionr = Append[Nsolutionr, {rx[nn], ry[nn]}];
  Nsolutionv = Append[Nsolutionv, {vx[nn], vy[nn]}];
  {vx[n], vy[n]} = {vx[nn], vy[nn]};
  {rx[n], ry[n]} = {rx[nn], ry[nn]}; , {n, 1, Floor[2500]}];
Dimensions[Nsolutionr];
ListPlot[{fullplot , Nsolutionr},  {PlotStyle -> {Blue, Red}, 
  PlotMarkers -> {〇, None},
  Joined -> {False, True}, PlotLabel -> "Newton Method", 
  ImageSize -> Large} ]

(*Second order Runge-Kutta*)

x0 = 10 ; v0 = 0.1 ; 
E0 = N[(1/2)*(v0)^2 - 1/x0 ]; a = N[-1/(2*E0)];
P = N[2*Pi*a^(3/2)] ; dt = P/500; 
rx[1] = 10; ry[1] = 0; vx[1] = 0; vy[1] = 0.1; 
RK2solutionr = {{rx[1], ry[1]}};
RK2solutionv = {{vx[1], vy[1]}}; n = 1;
Do[
  nn = n + 1;
  {ax[n], ay[n]} = - {rx[n], ry[n]}/(rx[n]^2 + ry[n]^2)^(3/2);
  {aax[n], 
    aay[n]} = - {rx[n] + dt*vx[n], 
      ry[n] + dt*
        vy[n]}/((rx[n] + dt*vx[n])^2 + (ry[n] + dt*vy[n])^2)^(3/2);
  {rx[nn], ry[nn]} = {rx[n], ry[n]} + {vx[n], vy[n]}*dt + 
    0.5*dt^2 * {ax[n], ay[n]};
  {vx[nn], vy[nn]} = {vx[n], vy[n]}  + 
    0.5*dt ({ax[n], ay[n]} + {aax[n], aay[n]})   ; 
  RK2solutionr = Append[RK2solutionr, {rx[nn], ry[nn]}];
  RK2solutionv = Append[RK2solutionv, {vx[nn], vy[nn]}];
  {rx[n], ry[n]} = {rx[nn], ry[nn]};
  {vx[n], vy[n]} = {vx[nn], vy[nn]};
  , {n, 1, 2500}] ;
Dimensions[RK2solutionr];
ListPlot[{fullplot, RK2solutionr},  {PlotStyle -> {Blue, Red}, 
  PlotMarkers -> {〇, None},
  Joined -> {False, True}, PlotLabel -> "2nd Order RK", 
  ImageSize -> Large} ]

(*Position Verlet Algorithm*)

x0 = 10 ; v0 = 0.1 ; 
E0 = N[(1/2)*(v0)^2 - 1/x0 ]; a = N[-1/(2*E0)];
P = N[2*Pi*a^(3/2)] ; dt = P/500; 
rx[1] = 10; ry[1] = 0; vx[1] = 0; vy[1] = 0.1; 
PVsolutionr = {{rx[1], ry[1]}};
PVsolutionv = {{vx[1], vy[1]}}; n = 1;
Do[
 nn = n + 1;
 {rx[nn], ry[nn]} = {rx[n], ry[n]} + 0.5*dt*{vx[n], vy[n]};
 {ax[n], ay[n]} = {rx[nn], ry[nn]}/(rx[nn]^2 + ry[nn]^2)^(3/2);
 {vx[nn], vy[nn]} = {vx[n], vy[n]} - dt*{ax[n], ay[n]};
 {rx[nn], ry[nn]} = {rx[nn], ry[nn]} + 0.5*dt*{vx[nn], vy[nn]}; 
 PVsolutionr = Append[PVsolutionr, {rx[nn], ry[nn]}];
 PVsolutionv = Append[PVsolutionv, {vx[nn], vy[nn]}];
 {vx[n], vy[n]} = {vx[nn], vy[nn]};
 {rx[n], ry[n]} = {rx[nn], ry[nn]};
 , {n, 1, 2500}]
Dimensions[PVsolutionr];
ListPlot[{fullplot, PVsolutionr},  {PlotStyle -> {Blue, Red}, 
  PlotMarkers -> {〇, None},
  Joined -> {False, True}, PlotLabel -> "Position Verlet", 
  ImageSize -> Large} ]


(*Plot energy error function*)

ce = Table[{n, (0.5*Norm[Csolutionv[[n]]]^2 - 
        1/Norm[Csolutionr[[n]]])/Abs[E0] + 1}, {n, 1, 2500 }];
ne = Table[{n, (0.5*Norm[Nsolutionv[[n]]]^2 - 
        1/Norm[Nsolutionr[[n]]])/Abs[E0] + 1}, {n, 1, 2500}];
rk2e = Table[{n, (0.5*Norm[RK2solutionv[[n]]]^2 - 
        1/Norm[RK2solutionr[[n]]])/Abs[E0] + 1}, {n, 1, 2500}];
pve = Table[{n, (0.5*Norm[PVsolutionv[[n]]]^2 - 
        1/Norm[PVsolutionr[[n]]])/Abs[E0] + 1}, {n, 1, 2500 }];

ListPlot[{Legended[ce, "Cromer"], Legended[ne, "Newton"], 
  Legended[rk2e, "RK2"], Legended[pve, "PV"]}, 
 PlotLabel -> "Energy Error Function Full Plot", ImageSize -> Large, 
 Joined -> True, Mesh -> Full, 
 Ticks -> {Table[{500*i, i}, {i, 1, 5}], Automatic}, 
 AxesLabel -> {"(t/P)", "(E/E0 + 1)"} ]

ce = Table[{n, (0.5*Norm[Csolutionv[[n]]]^2 - 
        1/Norm[Csolutionr[[n]]])/Abs[E0] + 1}, {n, 225, 275 }];
ne = Table[{n, (0.5*Norm[Nsolutionv[[n]]]^2 - 
        1/Norm[Nsolutionr[[n]]])/Abs[E0] + 1}, {n, 225, 275}];
rk2e = Table[{n, (0.5*Norm[RK2solutionv[[n]]]^2 - 
        1/Norm[RK2solutionr[[n]]])/Abs[E0] + 1}, {n, 225, 275}];
pve = Table[{n, (0.5*Norm[PVsolutionv[[n]]]^2 - 
        1/Norm[PVsolutionr[[n]]])/Abs[E0] + 1}, {n, 225, 275 }];
ListPlot[{Legended[ce, "Cromer"], Legended[ne, "Newton"], 
  Legended[rk2e, "RK2"], Legended[pve, "PV"]}, 
 PlotLabel -> "Energy Error Function at t/P range 0.45 to 0.55", 
 ImageSize -> Large, Joined -> True, Mesh -> Full, 
 Ticks -> {Table[{225 + i, N[0.45 + i/500]}, {i, 0, 50, 5}], 
   Automatic}, AxesLabel -> {"(t/P)", "(E/E0 + 1)"} ]

(*Find and plot angular momentum r cross v*)

AMC = Table[{n, 
    Csolutionr[[n, 1]]*Csolutionv[[n, 2]] - 
     Csolutionv[[n, 1]]*Csolutionr[[n, 2]]}, {n, 1, 2500}];

AMN = Table[{n, 
    Nsolutionr[[n, 1]]*Nsolutionv[[n, 2]] - 
     Nsolutionv[[n, 1]]*Nsolutionr[[n, 2]]}, {n, 1, 2500}];
AMRK2 = Table[{n, 
    RK2solutionr[[n, 1]]*RK2solutionv[[n, 2]] - 
     RK2solutionv[[n, 1]]*RK2solutionr[[n, 2]]}, {n, 1, 2500}];
AMPV = Table[{n, 
    PVsolutionr[[n, 1]]*PVsolutionv[[n, 2]] - 
     PVsolutionv[[n, 1]]*PVsolutionr[[n, 2]]}, {n, 1, 2500}];

ListPlot[{Labeled[AMC, {"Cromer Angular momentum"}], 
  Labeled[AMN, "Newton Angular Momentum"], 
  Labeled[AMRK2, "RK2 Angular Momentum"], 
  Labeled[AMPV, "Position Verlet Angular Momentum"]}, 
 PlotLayout -> "Row" , ImageSize -> Full, 
 Ticks -> {Table[{500*i, i}, {i, 1, 5}], Automatic}, 
 AxesLabel -> {"(t/P)", "(E/E0 + 1)"}]

ListPlot[{Legended[AMC, "Cromer Angular momentum"], 
  Legended[AMN, "Newton Angular momentum"], 
  Legended[AMRK2, "RK2 Angular momentum"], 
  Legended[AMPV, "PV Angular momentum"]}, 
 Joined -> True, {Ticks -> {Table[{500*i, i}, {i, 1, 5}], Automatic}, 
  AxesLabel -> {"(t/P)", "(E/E0 + 1)"}}]


(*Calculate v cross L - r*)

LC = Transpose[{ConstantArray[0, 2500], ConstantArray[0, 2500], 
    AMC[[All, 2]] }];
LN = Transpose[{{ConstantArray[0, 2500]}, {ConstantArray[0, 2500]}, 
    AMN[[All, 2]] }];
LRK = Transpose[{{ConstantArray[0, 2500]}, {ConstantArray[0, 2500]}, 
    AMRK2[[All, 2]] }];
LPV = Transpose[{{ConstantArray[0, 2500]}, {ConstantArray[0, 2500]}, 
    AMPV[[All, 2]] }];


cs3 = Table[{Csolutionr[[n, 1]], Csolutionr[[n, 2]], 
    ConstantArray[0, 2500]}, {n, 1, 2500} ];
ns3 = Table[{Csolutionr[[n, 1]], Nsolutionr[[n, 2]], 
    ConstantArray[0, 2500]}, {n, 1, 2500} ];
rk2s3 = Table[{Csolutionr[[n, 1]], RK2solutionr[[n, 2]], 
    ConstantArray[0, 2500]}, {n, 1, 2500} ];
pvs3 = Table[{Csolutionr[[n, 1]], PVsolutionr[[n, 2]], 
    ConstantArray[0, 2500]}, {n, 1, 2500} ];


AC = Cross[cs3, LC] - r; 
AN = Cross[ns3, LC] - r; 
ARK2 = Cross[rk2s3, LC] - r; 
APV = Cross[pvs3, LC] - r; 


