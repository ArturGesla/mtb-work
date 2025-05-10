postProcess -func writeCellCentres
grep rotor Ccx -A 4565   > a
grep rotor Ccy -A 4565   > b
paste a b > c 
python3 script.py
sed 's/^/(/g' vel | sed 's/$/)/' > vel2
cat p1 vel2 p3 > U
