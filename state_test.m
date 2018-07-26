
syms t real;
syms q1t(t) q2t(t);
syms q1 q2 q3 q4 dq1 dq2 dq3;

dq1t = diff(q1t,t);
dq2t = diff(q2t,t);

Test1 = [];
Test2 = [];

for q2 = 0:30
    for q3 = 0:30
        for q4 = 0:30
            Test1(q2+1,q3+1,q4+1) = -(2*(441326504143384259681916978611353695000*cos(2*q2) - 118587930864843491416670077701394965000*sin(2*q2) - 4885511350082437068281664924031155000*q3 - 7511816024975059696751642028140413000*cos(3*q2)*q3^2 + 3709792560642373878805177557456456000*sin(3*q2)*q3^2 - 587113327277922662700221950837093750*cos(q2)*q4 + 2686045782264848903731204906773562500*sin(q2)*q4 + 5278948389508569794772604158698411000*cos(q2)*q3^2 + 19396889719378550870424402320312148000*cos(q2)*q4^2 - 6505585656724111202406810344505828000*sin(q2)*q3^2 + 4239753669568830342970753826608494000*sin(q2)*q4^2 + 218782064667536329567054551718253059584*cos(2*q2)*q3*q4 + 273938936741414373271602956475160986624*sin(2*q2)*q3*q4 - 351601278404158657971154237063821165000))/(219573915450166586871606033693896247084*sin(2*q2) - 272214117513114385337688577191530939749*cos(2*q2) + 565394292314951515018344051549510676499);

            Test2(q2+1,q3+1,q4+1) = (7370471991418124945281175451205007066880*cos(3*q2) - 31104621833567842725626643674742867490560*sin(q2) - 13578659833993976454218532325997083226880*cos(q2) + 3284166492211939810514508194745155308800*sin(3*q2) + 78557584530055074470464879899053536000*q4 + 298764616417146195600645003683649185792*cos(2*q2)*q3^2 - 144270029255310751340024714461809441792*cos(4*q2)*q3^2 - 791850782630257304551481975643187500*cos(2*q2)*q4^2 + 374085789352246384550723455834330509312*sin(2*q2)*q3^2 - 32710060334843311687126287809268273920*sin(4*q2)*q3^2 + 1724819228299987933914379283630046875*sin(2*q2)*q4^2 + 199444139863663005811040890940078125*cos(q2)*q3 - 912457724579363794973452510912968750*sin(q2)*q3 - 37934562513766136418526238334451536000*cos(2*q2)*q4 + 30296539833823375403377676161749976000*sin(2*q2)*q4 + 21115793558034279179090416634793644000*cos(q2)*q3*q4 - 26022342626896444809627241378023312000*sin(q2)*q3*q4 - 30047264099900238787006568112561652000*cos(3*q2)*q3*q4 + 14839170242569495515220710229825824000*sin(3*q2)*q3*q4)/(219573915450166586871606033693896247084*sin(2*q2) - 272214117513114385337688577191530939749*cos(2*q2) + 565394292314951515018344051549510676499);
        end
    end
end