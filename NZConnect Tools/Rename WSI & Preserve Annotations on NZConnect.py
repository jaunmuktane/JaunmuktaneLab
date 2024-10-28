from getpass import getpass
from xml.etree import cElementTree as ET  # pip install elementpath
import requests
from dicttoxml import dicttoxml

# NZConnect parameters
# Login credentials
username = input("Username:\n")
password = getpass(prompt="Password:\n")

# Verify login credentials
nzconnectapi = 'https://nzconnect.ion.ucl.ac.uk/nz/connect/api/'
login = {'username': username, 'password': password}

try:
    login_request = requests.get(nzconnectapi + 'signin', params=login, timeout=60)
except requests.Timeout:
    print("NZConnect did not respond after 60 seconds. Request timed out.")

if 'ndperror' in login_request.text:
    print('Incorrect login credentials, please try again.')
    exit()
else:
    # Obtain sessionid for script access to NZConnect
    login_et = ET.fromstring(login_request.text)
    sessionid = login_et.find('sessionid').text


# Part 1 - Creating a copy of the image on the S drive
    #User selects the file to duplicate
    
    #User inputs the desired file name
    
    #Duplicating image on S drive


# Part 2 - Renaming duplicated image
    # Renaming duplicated file with file name input from user
    
    # Waiting until file appeis uploaded to NZconnect or running the nexr section separately


# Part 3 - Copying annotations on NZConnect

source_image_objectid = "D54C53F1-71E2-475F-801F-F820CE2399C9"
target_image_objectid = "1586A3B5-FCB7-4462-B4AE-349AAD0B2BC8"

source_anno_param = {'sessionid': sessionid, 'parentid': source_image_objectid, 'type': 'annotation'}
try:
    source_anno_request = requests.get(nzconnectapi + 'getchildren', params=source_anno_param, timeout=60)
except requests.Timeout:
    print("NZConnect did not respond after 60 seconds. Request timed out.")

# Check if image has annotations
if 'resultsreturned="0"' in source_anno_request.text:
    print('No annotations found for this image.')
    exit()
else:
    source_anno_et = ET.fromstring(source_anno_request.text)

# Collect all annotation objectids
source_anno_objectids = [anno.find('id').text for anno in source_anno_et.findall('object')]

# Get annotation XML (contain annotation coordinates)
anno_xml = []
for anno_objectid in source_anno_objectids:
    anno_object_param = {'sessionid': sessionid, 'objectid': anno_objectid}

    try:
        anno_object_request = requests.get(nzconnectapi + 'getobject', params=anno_object_param, timeout=60)
    except requests.Timeout:
        print("NZConnect did not respond after 60 seconds. Request timed out.")
    
    anno_xml.append(anno_object_request.text)

# Add annotations to target image individually using POST request
for xml in anno_xml:
    # Annotation parameters
    name = ET.fromstring(xml).find('name').text
    type = ET.fromstring(xml).find('type').text
    inherit = True
    dateadded = ET.fromstring(xml).find('dateadded').text
    owner = ET.fromstring(xml).find('owner').find('id').text
    content = ET.tostring(ET.fromstring(xml).find('content'), encoding='unicode', method='xml')
    
    add_anno_param = {
        'sessionid': sessionid,
        'parentid': target_image_objectid,
        'name': name,
        'type': type,
        'inherit': inherit,
        'dateadded': dateadded,
        'owner': owner,
        'content': content
    }

    # Convert the parameter dictionary to an XML string
    xml_data = dicttoxml(add_anno_param, custom_root='annotation', attr_type=False).decode('utf-8')

    # Send POST request to add annotation
    try:
        add_anno_request = requests.post(
            nzconnectapi + 'addobject',
            data=xml_data,
            headers={'Content-Type': 'text/xml'},
            timeout=60
        )
        add_anno_request.raise_for_status()
    except requests.Timeout:
        print("NZConnect did not respond after 60 seconds. Request timed out.")
