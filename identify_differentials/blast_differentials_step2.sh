
cat "$@" | sort -g | uniq -c > tmp.csv

while read line; do
#for LINE in 'tmp.csv'; do
#	echo "$LINE"
	echo "$line" | grep 'differentials.csv' 
	COPIES=$(echo "$line" | cut -d' ' -f1)
	ID=$(echo "$line" | cut -d' ' -f2)
	if [ "$COPIES" != "$#" ]; then
#		echo "copies are: $COPIES"
#		echo "ID is: $ID"
#		echo		
		printf "$ID" >> differentials.csv
		printf "\n" >> differentials.csv
		grep "$ID" differentials_tab.csv 
	else
		printf "$ID" >> in_all.csv
		printf "\n" >> in_all.csv		
	fi
done<'tmp.csv'
 
rm tmp.csv 

exit