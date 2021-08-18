import subprocess
import re

# checkm_taxon_file obtained with command 'checkm taxon_list' in checkM env
checkm_taxon_file='files/checkm_taxon_list.txt'
checkm_genera_table = subprocess.check_output(['grep', 'genus', checkm_taxon_file])
checkm_genera = re.findall(r'\w+', str(checkm_genera_table))

genus_name=False
for word in checkm_genera:
    if genus_name:
        print(word)
    if word == 'genus':
        genus_name=True
    else:
        genus_name=False