
SETS
n /n1*n6/
g /g1*g5/
d /d1*d4/
l /l1*l9/
pros(l) /l4*l9/
ex(l) /l1*l3/
mapG(g,n) /g1.n1 ,g2.n2 ,g3.n3 ,g4.n5 ,g5.n6/
mapD(d,n) /d1.n1 ,d2.n4 ,d3.n5 ,d4.n6/
ref(n) /n1/
mapSL(l,n) /l1.n1 ,l2.n1 ,l3.n4 ,l4.n2 ,l5.n2 ,l6.n3,l7.n3,l8.n4,l9.n5/
mapRL(l,n) /l1.n2 ,l2.n3 ,l3.n5 ,l4.n3 ,l5.n4 ,l6.n4,l7.n6,l8.n6,l9.n6/

blake1 /'B', 'FLmax', 'IC'/
;

TABLE LDATA(l, blake1)
    B   FLmax    IC
l1  500 150      0
l2  500 150      0
l3  500 150      0
l4  500 150      700000
l5  500 200      1400000
l6  500 200      1800000
l7  500 200      1600000
l8  500 150      800000
l9  500 150      700000;

SCALAR ILmax
/3000000/;

TABLE DDATA(d,*)
         PDmax   LScost
d1       200     40
d2       150     52
d3       100     55
d4       200     65;
TABLE GDATA(g,*)
        PEmax   Gcost
g1      300     18
g2      250     25
g3      400     16
g4      300     32
g5      150     35;

SCALAR SIGMA
/8760/;
SCALAR M
/5000/;
VARIABLES
Z
PL(l)
THETA(n);
POSITIVE VARIABLES
PG(g)
PLS(d);
BINARY VARIABLES
x(l);

EQUATIONS
EQ3A
EQ3B
EQ3D
EQ3E
EQ3Fa
EQ3Fb
EQ3Ga
EQ3Gb
EQ3Ha
EQ3Hb
EQ3I
EQ3J
EQ3Ka
EQ3Kb
EQ3L
;

EQ3A.. Z=E=SUM(l$pros(l),LDATA(l,'IC')*x(l))+SIGMA*(SUM(g,GDATA(g,'Gcost')*PG(G))+SUM(d,DDATA(d,'LSCOST')*PLS(d)));
EQ3B.. SUM(l$pros(l),LDATA(l,'IC')*x(l)) =L= ILmax;

*EQUATIONS 3C ARE DEFINITION OF BINARY VARIABLES

EQ3D(n).. SUM(g$mapG(g,n),PG(G))-SUM(l$mapSL(l,n),PL(l))+SUM(l$mapRL(l,n),PL(l))=E= SUM(d$mapD(d,n),DDATA(d,'PDMAX')-PLS(d));
EQ3E(l)$EX(l).. PL(l)=E=LDATA(l,'B')*(SUM(n$mapSL(l,n),THETA(n))-SUM(n$mapRL(l,n),THETA(n)));
EQ3Fa(l)$EX(l).. -LDATA(l,'FLmax')=L=PL(l);
EQ3Fb(l)$EX(l).. PL(l)=L=LDATA(l,'FLmax');
EQ3Ga(l)$PROS(l).. -x(l)*LDATA(l,'FLmax')=L=PL(l);
EQ3Gb(l)$PROS(l).. PL(l)=L=x(l)*LDATA(l,'FLmax');
EQ3Ha(l)$PROS(l).. -(1-x(l))*M=L=PL(l)-LDATA(l,'B')*(SUM(n$mapSL(l,n),THETA(n))-SUM(n$mapRL(l,n),THETA(n)));
EQ3Hb(l)$PROS(l).. PL(l)-LDATA(l,'B')*(SUM(n$mapSL(l,n),THETA(n))-SUM(N$mapRL(l,n),THETA(n)))=l=(1-X(l))*M;
EQ3I(g).. PG(g)=L=GDATA(g,'PEmax');
EQ3J(d).. PLS(d)=L=DDATA(d,'PDmax');
EQ3Ka(n).. -3.14=L=THETA(n);
EQ3Kb(n).. THETA(n)=L=3.14;
EQ3L(n)$REF(n).. THETA(n)=L=0;

MODEL TEP_DET /ALL/;
SOLVE TEP_DET USING MIP MINIMIZING Z;
