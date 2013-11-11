#!/bin/bash
rsync -avz cz1@129.109.88.204:lwork/vir/ZrD2r1n20.dat .
rsync -avz cz1@129.109.88.204:lwork/vir/ZrD3r1n16.dat .
rsync -avz cz1@129.109.88.204:lwork/vir/ZrD4r1n12.dat .
rsync -avz cz1@129.109.88.204:lwork/vir/ZrD5r1n16.dat .
rsync -avz cz1@129.109.88.204:lwork/vir/ZrD6r1n20.dat .
rsync -avz cz1@129.109.88.204:lwork/vir/ZrD7r1n24.dat .
rsync -avz cz1@129.109.88.204:lwork/vir/ZrD8r1n28.dat .
python mkZ.py
