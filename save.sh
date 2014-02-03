#read -p "Enter prob value:" i
cp *.his kprob0/
cp HISTG* kprob0/
cp screen kprob0/
cp fort.* kprob0/
cp -rf design/ kprob0/.
cp -rf results/ kprob0/.
echo "Success saving the results into kprob"
