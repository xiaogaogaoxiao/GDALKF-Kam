clear all;
close all;
clc;

syms ax ay az real
syms mx my mz real
syms mN mD real
syms q0 q1 q2 q3 real
syms mu real

q = [q0; q1; q2; q3];

C = q2c(q);

Ab = [ax; ay; az];
Mb = [mx; my; mz];


Ar = [0; 0; 1];
Mr = [mN; 0; mD];


fA = Ab - C * Ar;

fAM = [Ab; Mb]-[C * Ar; C * Mr];



JA = jacobian(fA, q);

JAM = jacobian(fAM, q);



delta_qA = simplify( - JA' * fA);

delta_qAM = simplify( - JAM' * fAM);



qA = q + mu * delta_qA;

qAM_ = q + mu * delta_qAM;

qAM = [
 q0 - mu*(4*q0 - 2*mz*mD*q0 - 2*my*mD*q1 + 2*mx*mD*q2 - 2*mx*mN*q0 - 2*mz*mN*q2 + 2*my*mN*q3 - 2*az*q0 - 2*ay*q1 + 2*ax*q2);
 q1 - mu*(4*q1 - 2*my*mD*q0 + 2*mz*mD*q1 - 2*mx*mD*q3 - 2*mx*mN*q1 - 2*my*mN*q2 - 2*mz*mN*q3 - 2*ay*q0 + 2*az*q1 - 2*ax*q3);
 q2 - mu*(4*q2 + 2*mx*mD*q0 + 2*mz*mD*q2 - 2*my*mD*q3 - 2*mz*mN*q0 - 2*my*mN*q1 + 2*mx*mN*q2 + 2*ax*q0 + 2*az*q2 - 2*ay*q3);
 q3 - mu*(4*q3 - 2*mx*mD*q1 - 2*my*mD*q2 - 2*mz*mD*q3 + 2*my*mN*q0 - 2*mz*mN*q1 + 2*mx*mN*q3 - 2*ax*q1 - 2*ay*q2 - 2*az*q3)
];

Kam = jacobian(qAM, q)

simplify(qAM - Kam * q)