 cat run1 | grep "for Ux" | cut -f 2 -d "=" | cut -f 1 -d ',' > resUx
 cat run1 | grep "for Uy" | cut -f 2 -d "=" | cut -f 1 -d ',' > resUy
 cat run1 | grep "for Uz" | cut -f 2 -d "=" | cut -f 1 -d ',' > resUz
 cat run1 | grep "for p" | cut -f 2 -d "=" | cut -f 1 -d ',' > resP
# cat run1 | grep "for k" | cut -f 2 -d "=" | cut -f 1 -d ',' > resK
# cat run1 | grep "for eps" | cut -f 2 -d "=" | cut -f 1 -d ',' > resEps

#gnuplot -e 'set grid; set logscale y; plot "resUx" w l,"resUy" w l,"resUz" w l,"resP" w l,"resK" w l ,"resEps" w l; pause -1'
gnuplot -e 'set grid; set logscale y; plot "resUx" w l,"resUy" w l,"resUz" w l,"resP" w l; pause -1'

