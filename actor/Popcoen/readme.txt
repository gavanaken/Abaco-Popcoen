************************************************
* Documentation how to install and run Popcoen *
************************************************

1. Popcoen Installation

As you read this text, you have extracted all popcoen files to one folder and opened this readme.txt file. 

Here, it is explained how to get Popcoen running. 
We assume that Linux, Python >= 2.7 (not Python 3.x) and pip is installed.

* go to the folder of this readme.txt file by
cd [folder]

* make sure that all .py files are executable, by
chmod u+x *.py

* install mdtraj with 
pip install mdtraj --upgrade

* install the Popcoen neural network with
pip install PopcoenNN_package.zip

* make sure that you have the actual versions of the following packages 
numpy, scipy, scikit-learn, theano. Just try updating them with
pip install [package] --upgrade

Installation completed. 

-------------------------------------------------------------------------------------------

2. Popcoen application

2a. Running Popcoen without server

The easiest way to make an entropy prediction is by piping a protein structure (in pdb format) into popcoenMain.py, for example like this:

cat 3x2m.pdb | ./popcoenMain.py

This will load Popcoen into memory, process the structure, make the entropy prediction, and finish Popcoen.
As loading Popcoen takes several seconds, this way is not suitable for many structures. 
However, it can be useful for very limited number of structures (avoiding server-client communication). 

- - - - - - - - - - - - - - - - - - - - - - -

2b. Running Popcoen as server

For entropy prediction  of many structures, run Popcoen as a server. 
Then, Popcoen is loaded to memory only once when the server is started. 
Proceed as follows:

* start the Popcoen server via
./popcoen_server.py

This starts a server which listens on the 10888 port for incoming messages. 
If this port is already in use, please change the server_address line of popcoen_server.py.
By default, the server is reachable on localhost. 
This can easily be changed to the Internet (see popcoen_server.py).

* the server waits for incoming requests. This requests can come from any software, 
which allows to include Popcoen easily in existing software just by establishing socket communication. 
We also provide a client for the server. Open another terminal and call:
cat 3x2m.pdb | ./client_for_popcoen.py

This should give the same prediction as in 2a, however achieved much quicker.

--------------------------------------------------------------------


Thank you for your interest in Popcoen. 

For any questions, please contact Martin Goethe via martingoethe@gmx.de. 
Also, if you try to incorporate Popcoen into your software, please let us know.
There might be a tailored way to improve prediction precision for your application. 

If you use Popcoen, please cite our article describing Popcoen. 

Best regards,
Martin Goethe    (Barcelona, May 1st, 2017)
