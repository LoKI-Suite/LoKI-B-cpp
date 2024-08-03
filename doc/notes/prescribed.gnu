g(s,r)=sqrt(r)*s*gamma(2.5/s)**1.5*gamma(1.5/s)**-2.5*exp(-(gamma(2.5/s)/gamma(1.5/s)*r)**s)
set xlabel "r"
set ylabel "g_s(r)"
set term pdf color
plot [r=0:3] g(0.5,r) title "s=0.5", g(1,r) title "s=1", g(2,r) title "s=2", \
	g(5,r) title "s=5", g(10,r) title "s=10", g(20,r) title "s=20", \
	g(50,r) title "s=50", g(100,r) title "s=100", g(1000,r) title "s=1000"
