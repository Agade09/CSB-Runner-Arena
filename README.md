# CSB-Runner-Arena
A local arena for Coders Strike Back (codingame.com) which tries to gauge the performance of runners in a collisionless environment. Only works on Linux.

This program will test the speed of your runners by having it run 100 maps taken from the CG IDE. At the end it outputs the total number of turns your AI took to finish all 100 races without an opponent in a collisionless simulation (your runners don't collide with each other).

Usage:
* Compile the Arena using the Makefile
* Compile your AI (as is from CodinGame)
* Run the Arena with the name of your AI's binary as a parameter (e.g: Arena MyAI-V23)
