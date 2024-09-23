# Environment setup for Cell Ranger

This section is primarily based on the instructions provided by this website: [setup-and-install-cellranger](https://labs.wsu.edu/winuthayanon/how-to-analyze-single%E2%80%90cell-rna%E2%80%90seq/setup-and-install-cellranger/). 


1. Download and unpack the Cell Ranger file in any directory. In this example, we will use `$HOME/opt`:

	```
	mkdir -p $HOME/opt
	cd $HOME/opt
	```

2. Download cellranger from [10xgenomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads):

	```
	wget -O cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1725683732&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=eOPXD1vEBDioK~mLi39txS90verKWBrbq-45dr9mZ61t1gbJXW4q3DIjRy58EXxs4yI6E8~xbvePYd6OzvALcsXyUGQfDZkWbbxCM~J9BlIyZxCJcdkgdzBT1w7uTHpar3WnwfgsHT70VVVMjRBu7r8MglmVW7oHUqZGEVsSqRma3J5GSjN6uwKmfsXkyds4KHgAhFIuFnZGKLUH6~akARXPpy2iKKIiHyV8vui5~Ond0VSq1QQtJMRwnbO4CzPPT~Wmg5Qi9wI4-UmN~d1QZkM~VElrAj~Va3GcNitL6~boalAXnTKTUBoUM5UpEPI-u5mY0DA5jRumydfN1y6xFQ__"
	```

	Then, unpack the downloaded file:

	```
	tar -xzvf cellranger-8.0.1.tar.gz
	```

3. Before running the Cell Ranger commands, make sure to prepend the Cell Ranger directory to your `$PATH`, allowing you to invoke the cellranger command from any location: 

	```
	export PATH=$HOME/opt/cellranger-8.0.1:$PATH
	```

	***You may want to add this command to your `$HOME/.bashrc` for convenience.***