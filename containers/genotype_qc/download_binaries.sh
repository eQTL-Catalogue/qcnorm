#Make bin directory
mkdir bin

#Download GenotypeHarmonizer
wget http://www.molgenis.org/downloads/GenotypeHarmonizer/GenotypeHarmonizer-1.4.20-dist.tar.gz
tar xzfv GenotypeHarmonizer-1.4.20-dist.tar.gz
mv GenotypeHarmonizer-1.4.20-SNAPSHOT/ bin/GenotypeHarmonizer-1.4.20
rm GenotypeHarmonizer-1.4.20-dist.tar.gz

#Download LDAK
wget http://dougspeed.com/wp-content/uploads/ldak5.linux_.zip
unzip ldak5.linux_.zip
mv ldak5.linux bin/ldak
rm ldak5.linux_.zip
