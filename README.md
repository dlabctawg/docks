# Welcome to the Docks!

![Is that blue seal a chipmunk?](https://blog.docker.com/media/2015/09/animals-august2015.png)

The Docks at the D-Lab Computational Text Analysis Working Group is a series of projects based off of a common docker image. Try it yourself! With Docker installed and the docker deamon running, try this in your command line:

```
cd ~ # Or navigate to your favorite GitHub project directory!
git clone https://github.com/dlabctawg/docks
cd docks
docker run -v $(pwd):/home/rstudio --name rstudio -d -p 8080:8787 dlabctawg/docks
open http://localhost:8080/ # or open that URL in a browser window
```

An Rstudio Server login window should pop up. Type `rstudio` for the user and `rstudio` for the password, and you should see the project directories under files. Open a directory and start playing with project resources!