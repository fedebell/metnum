#!/bin/bash
#Connesione a phoenix
ssh fbelliardo@192.168.80.24
expect "fbelliardo@192.168.80.24's password: "
send "$pass"
screen -dm -S 2313 ./sg -p 0.2313 4 6 8 10 12 14 16 18 20 22 24
screen -dm -S 2361 ./sg -p 0.2361 4 6 8 10 12 14 16 18 20 
screen -dm -S 2340 ./sg -p 0.2340 4 6 8 10 12 14 16 18 20 22
exit
#Connesione a defiant
ssh fbelliardo@192.168.80.27
expect "fbelliardo@192.168.80.27's password: "
send "$pass"
screen -dm -S 2410 ./sg -p 0.2410 4 6 8 10 12 14 16 18
screen -dm -S 2420 ./sg -p 0.2420 4 6 8 10 12 14 16 
screen -dm -S 2440 ./sg -p 0.2440 4 6 8 10 12 14 16 
exit
#Connesione a reliant
ssh fbelliardo@192.168.80.31
expect "fbelliardo@192.168.80.31's password: "
send "$pass"
screen -dm -S 2460 ./sg -p 0.2460 4 6 8 10 12 14
screen -dm -S 2480 ./sg -p 0.2480 4 6 8 10 12 14 
screen -dm -S 2500 ./sg -p 0.2500 4 6 8 10 12 14
exit
#Connessione a saratoga
ssh fbelliardo@192.168.80.33
expect "fbelliardo@192.168.80.33's password: "
send "$pass"
screen -dm -S 2520 ./sg -p 0.2520 4 6 8 10 12 14
screen -dm -S 2540 ./sg -p 0.2540 4 6 8 10 12 14
screen -dm -S 2560 ./sg -p 0.2560 4 6 8 10 12 14
exit
#Connessione a intrepid
ssh fbelliardo@192.168.80.24
expect "fbelliardo@192.168.80.24's password: "
send "$pass"
screen -dm -S 2580 ./sg -p 0.2580 4 6 8 10 12 14
screen -dm -S 2600 ./sg -p 0.2600 4 6 8 10 12 
screen -dm -S 2620 ./sg -p 0.2620 4 6 8 10 12 
exit
#Connessione a coccoina
#ssh fbelliardo@192.168.111.91
#expect "fbelliardo@192.168.111.91's password: "
#send "$pass"
#screen -dm -S 2248 ./sg -p 0.2248 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32
exit
 
