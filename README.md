Designing primers for regions surrounding guide rna

## Usage

To install on dockerfile ensure docker and docker compose is installed on the system

Run Docker Compose:
    ```
    docker compose up
    ```
    The jupyter notebook server should be running at `localhost:8888`. and the terminal should have the key to access the server. This docker enviornment ensures that all requirments are there to run the software.


Alternatively you can copy `create_amplicons.py` into an already working environment

to run the software 
```
python create_amplicons.py guides.tsv
```

and the output should be in output.csv
