#read -p "Enter prob value:" i
mv *.his kprob0/
mv HISTG* kprob0/
mv screen kprob0/
mv fort.* kprob0/
cp iterate.dat kprob0/.
cp -rf design/ kprob0/.
cp -rf results/ kprob0/.
echo "Success saving the results into kprob0"
