from PhaLP import Logger, get_url, download_VirusHostDB, download_NCBI_taxonomy, download_Expasy_enzyme, read_json, save_json, update
import sys
from datetime import datetime
import os
import feedparser
from pathlib import Path
from time import sleep

sys.stdout = Logger()
url_html = get_url()
url_txt = url_html+'&format=txt'
# change directory to map above 'code', since all other files are in there
os.chdir('/'.join(os.getcwd().split('/')[:-1]))

if Path('picklejar/last_updated.json').exists():
    last_PhaLP_update, last_PhaLP_version = read_json('picklejar/last_updated.json')
    last_PhaLP_update = datetime.strptime(last_PhaLP_update,'%a, %d %b %Y %H:%M:%S %Z')
else:
    last_PhaLP_update = datetime.strptime('13 Nov 2019 00:00:00 GMT','%d %b %Y %H:%M:%S %Z')
    last_PhaLP_version = '2019_10'
print('Last PhaLP version: %s (%s)'%(last_PhaLP_version,last_PhaLP_update))

while True:
    try:
        d = feedparser.parse(url_txt, modified=last_PhaLP_update)
        if d.status == 200:
            PhaLP_version = d['headers']['X-UniProt-Release']
            print('New PhaLP version: %s (%s)'%(PhaLP_version,datetime.strptime(d.modified,'%a, %d %b %Y %H:%M:%S %Z')))

            break
    except AttributeError as err:
        print('\n',err, ': retry')


print('#UPDATE downloading files ')
download_VirusHostDB()
download_NCBI_taxonomy()
download_Expasy_enzyme()               
print('#UPDATE analyzing PhaLPs ')

analyzed = update(last_PhaLP_version, PhaLP_version, url_html)
a = read_json('picklejar/analyzed.json')   
if last_PhaLP_version in a.keys():
    old_phalps = a[last_PhaLP_version]['phalps']
else: 
    old_phalps = []
new_phalps = []
for acc in analyzed['phalps']:
    if acc not in old_phalps:
        new_phalps.append(acc)

last_PhaLP_update = d.modified
last_PhaLP_version = PhaLP_version
save_json([last_PhaLP_update, last_PhaLP_version],'picklejar/last_updated.json')
print('%s new Phage Lytic Proteins are added in the latest PhaLP update.(%s PhaLPs total)'%(len(new_phalps),len(analyzed['phalps'])))
