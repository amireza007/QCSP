$title A Model for scheduling QCs
$onEolCom

scalars
    a1 "Weight for the makespan" /1/
    a2 "Weight for the total completion time" /1/
    M /1e+14/
;

Sets
    i "index of tasks" /0*20/
    dummy "final position index" /'0',T/
    k "index of QCs" /1*4/
    
    psi(i,i) "set of pairs (i,j) that can't performed simultaneously" / 1.(2*15),
                                                                        2.(3*15),
                                                                        3.(4*15)
                                                                        4.(5*15),
                                                                        5.(6*15),
                                                                        6.(7*15),
                                                                        7.(8*15),
                                                                        8.(9*15),
                                                                        9.(10*15),
                                                                        10.(11*15),
                                                                        11.(12*15)
                                                                        12.(14*15),
                                                                        13.(14*15),
                                                                        14.15/
    phi(i,i) "precedence relation of two tasks" /12.13/
    YS(*) "set of QC locations" /1*20/

;

alias(i,j);
alias(v,k);
alias(u,i);

Parameters
    p(i)
    r(k)
    l(i) "location of task i" / 1 1,
                                2 2,
                                3 4,
                                4 4,
                                5 6,
                                6 11,
                                7 12,
                                8 13,
                                9 14,
                                10 15,
                                11 16,
                                12 18,
                                13 18,
                                14 19,
                                15 20 /
    
    l0(k) "starting location of QC k" /1 1, 2 3, 3 5, 4 7/
    lT(k) "final position of QC k" /1 1, 2 3, 3 5, 4 7/ !!I just assusmed that QCs return to their initial position after their last jobs
    t(i,j) "Travel time of a QC from locationn li of task i to location lj of task j."
    t_d(dummy,i,k)
;

t(i,j) = abs(l(i)-l(j));
p(i) = uniform(60,90);
r(k) = 0;
t_d('T',i,k) =abs(l(i) - lt(k));
t_d('0',i,k) = abs(l(i) - l0(k));

Binary Variables
    X(i,j,k)    "1, if QC k performs task j IMMeDIATELY after performing i"
    Z(i,j)      "1, if task j starts later than the completion time of task i;"
    X_d(dummy, i ,k)
;

Positive Variables
    Y(k) "completion time of QC k"
    D(i) "completion time of task i"
    W "time at which all tasks are completed"
;
variable objFunc;

Equations
    obj, c1(k), c2(k)
    ,c3(k)
    ,c4(j)
    ,c5(i,k)
    ,c6(i,j,k)
    ,c7(i,j)
    ,c8(i,j)
    ,c9(i,j)
    ,c10(i,j,k)
    ,c11(j,k)
    ,c12(j,k)
;

obj.. objFunc =e= a1*W + a2*sum(k,Y(k));

c1(k).. Y(k) =l= W;

c2(k).. sum(j,x_d('0',j,k)) =e=1;

c3(k).. sum(i,x_d('T',i,k)) =e= 1;

c4(j).. sum((k,i)$(j.val <> i.val), x(i,j,k))  =e= 1;

c5(i,k).. sum(j $(j.val <> i.val),x(i,j,k)) - sum(j$(j.val <> i.val),x(j,i,k)) =e= 0;

c6(i,j,k).. D(i) + t(i,j) + p(j) - D(j) $(j.val <> i.val) =l= M*(1- x(i,j,k));

c7(i,j)..D(i) + p(j) $(phi(i,j) and (j.val <> i.val))  =l= D(j);

c8(i,j).. D(i) - D(j) + p(j)$(j.val <> i.val)  =l= M*(1 - z(i,j));

c9(i,j).. Z(i,j) + z(j,i) $(psi(i,j)and (j.val <> i.val)) =e= 1;

c10(i,j,k).. sum((v,u)$(v.val <= k.val), x(u,j,v)) -sum((v,u)$(v.val <= k.val), x(u,i,v)) $(j.val <> i.val and l(i)<l(j))  =l= M*(z(i,j) + z(j,i));

c11(j,k).. D(j) + t_d('T',j,k) - Y(k) =l= M*(1 - x_d('T',j,k));

c12(j,k).. r(k) - d(j) + t_d('0',j,k) + p(j) =l= M*(1 - x_d('0',j,k));


Model QCSP /all/;
solve QCSP using MINLP min objFunc;







