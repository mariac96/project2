ΚΑΡΛΗΣ ΝΙΚΟΛΑΟΣ 1115201800068
ΧΑΤΖΗΠΑΥΛΟΥ ΜΑΡΙΑ 1115201400223

Το Ai εκτελείται με τις εντολές όταν δίνουμε για algorithm  LSH ή Hypercube π.χ. με:
./search -i nasdaq2015_2017.csv -q .csv -M 10 -probes 2 -o outputfile2 -algorithm Hypercube

./search -i nasdaq2015_2017.csv -q .csv -o output_file -algorithm LSH


Το Ai είναι το LSH ή Hypercube της προηγούμενης εργασίας με ελάχιστες αλλαγές.

Στο Aii δίνεται algorithm Frechet και metric discrete. Στην συνάρτηση discreteFrechet μετατρέπουμε το dataset από την μορφή (Y1, Y2,..., Yn) σε μορφή [(1,Υ1), (2,Υ2), ..., (n, Υn)] και τα αποθηκεύουμε με μια κλάση curve που περιέχει ένα vector με points(x, y).
Σε κάθε curve εφαρμόζουμε το grid, και στην νέα καμπύλη αφαιρούμε τα duplicates και κάνουμε concat σε ενα 1d vector όπου και κάνουμε το padding. Eφαρμόζουμε την συνάρτηση g(p) της προηγούμενης εργασίας και βρίσκουμε το bucket. Επαναλαμβάνουμε L φορές.

Στο Aiii δίνεται algorithm Frechet και metric continuous. Στην συνάρτηση discreteontinuous μετατρέπουμε το dataset από την μορφή (Y1, Y2,..., Yn) σε μορφή [(1,Υ1), (2,Υ2), ..., (n, Υn)] και τα αποθηκεύουμε με μια κλάση curve που περιέχει ένα vector με points(x, y).
Μετατρέπουμε και τα δεδομένα που στις καταλληλες κλάσεις της βιβλιοθήκης fred.
Σε κάθε curve κάνουμε το filtering.
Μετά εφαρμόζουμε το snapping στο grid, μετά minima maxima, και κάνουμε το padding. Κάνουμε hash και μετά αποθηκεύουμε στο κατάλληλο bucket.
