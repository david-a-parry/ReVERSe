# ReVERSe

**R**ar**e** **V**ariant **E**nrichment and **R**ecessive **Se**gregation

### Installation
    
    git clone https://git.ecdf.ed.ac.uk/dparry/reverse.git
    cd reverse
    python3 -m pip install . --process-dependency-links

If you get an error "no such option: --process-dependency-links" (with newer
versions of pip) with the third command, the following should work:

    python3 -m pip install -r requirements.txt
    python3 -m pip install .

Remember to add the --user flag to pip commands if you get permissions errors.


