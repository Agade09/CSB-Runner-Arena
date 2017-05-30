# CSB-Runner-Arena
A local arena for Coders Strike Back (codingame.com) which tries to gauge the performance of runners in a collisionless environment. **Linux Only**.

This program will test the speed of your runners by having it run against maps generated the same way Codingame generates them. Arena outputs average turns taken to complete maps.

## Usage:
* Compile the Arena using the Makefile
* Compile your AI
* Run the Arena with the name of your AI's binary as a parameter (e.g: Arena MyAI-V23)

## Optional:
* Specify the number of threads as a command line parameter. e.g: Arena MyAI-V23 2
* Set timeout behavior on or off via the "constexpr bool Timeout" variable. This can be useful as I've noticed timeouts if the computer is being used for something else.
