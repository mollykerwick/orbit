# orbit
A simple model of the orbits of Earth and Jupiter around the sun. Model in C++ with a simple plot in Python. You can use this model to simulate what would happen to the orbit of these planets if the mass of the planets and/or the mass of the sun changed.

Needed: Pandas and Matplotlib.pyplot



To run this simple simulation:

`g++ -g orbit.cpp -o orbit -Wall`

`./orbit`

`python3 plot.py`



=====MINICONDA INSTALLATION ISSUES=====

from: https://stackoverflow.com/questions/31615322/zsh-conda-pip-installs-command-not-found

SUMMARY:

For the miniconda installation to initialize upon opening the terminal, add `source ~/.bash_profile` to your ~/.zshrc. This is because the miniconda installation instructions call for editing your PATH from ~/.bash_profile, but for macOS, the terminal default runs zsh instead of bash.
