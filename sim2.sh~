#!/bin/bash
#Connesione a phoenix
set pass "feupaecaim"
ssh fbelliardo@192.168.80.24
expect "fbelliardo@192.168.80.24's password: "
send "$pass"
screen -dm -S 224 ./sg -p 0.224 34 36 38 40
screen -dm -S 2265 ./sg -p 0.2265 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 
exit
#Connesione a defiant
ssh fbelliardo@192.168.80.27
expect "fbelliardo@192.168.80.27's password: "
send "$pass"
screen -dm -S 2391 ./sg -p 0.2391 16 18 20 22 24 26
screen -dm -S 2293 ./sg -p 0.2293 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34
exit
#Connesione a reliant
ssh fbelliardo@192.168.80.31
expect "fbelliardo@192.168.80.31's password: "
send "$pass"
screen -dm -S 2327 ./sg -p 0.2327 22 24 26 28 30 32
screen -dm -S 2300 ./sg -p 0.2300 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32
exit
#Connessione a saratoga
ssh fbelliardo@192.168.80.33
expect "fbelliardo@192.168.80.33's password: "
send "$pass"
screen -dm -S 2275 ./sg -p 0.2275 32 34 36
screen -dm -S 2350 ./sg -p 0.2350 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 
exit
#Connessione a intrepid
ssh fbelliardo@192.168.80.24
expect "fbelliardo@192.168.80.24's password: "
send "$pass"
screen -dm -S 2255 ./sg -p 0.2255 34 36 38 40
screen -dm -S 2370 ./sg -p 0.2370 4 6exit 8 10 12 14 16 18 20 22 24 26 28 30 32
exit
#Connessione a coccoina
#ssh fbelliardo@192.168.111.91
#expect "fbelliardo@192.168.111.91's password: "
#send "$pass"
screen -dm -S 2248 ./sg -p 0.2248 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32
exit
 
