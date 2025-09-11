import sys
import xml.etree.ElementTree as ET

def replace_testsuite_name(filename, new_name, target_index):
    tree = ET.parse(filename)
    root = tree.getroot()

    match_index = 0
    for child in root:
        if match_index == target_index:
            if child.tag == "testsuite" and child.get('name', "").strip() == "testdrive":  
                child.set("name", new_name)
                break
        match_index += 1

    for suite in root.findall("testsuite"):
        if suite.get("name") == new_name:
            for testcase in suite.findall("testcase"):
                classname = testcase.get("classname")
                testcase.set("classname", new_name)

    tree.write(filename, encoding="utf-8", xml_declaration=True)
    
if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: python3 script.py <xmlfile> <new_name> <index>")
    filename = sys.argv[1]
    new_name = sys.argv[2]
    try:
        target_index = int(sys.argv[3])
    except ValueError:
        sys.exit("Index must be an integer")
    replace_testsuite_name(filename, new_name, target_index - 1)
