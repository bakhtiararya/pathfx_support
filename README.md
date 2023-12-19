# PathFX Spport

Author: Arya Bakhtiar


## Google Drive Folder Link

Due to memory constraints on Github, all Inputs and Output Folders are posted to the following Google Drive folder. Must require permission to access the folder.

## About The Project

This project was meant to provide analytical and visual data regarding the STRING and PathFX protein interaction networks. The information provided by this repo is presented in graphs (.png), text (.txt), and tables (.csv) files regarding the following:

* Shared edges (shared protein/drug interactions present) 
* Shared Nodes (shared proteins/drug present)
* Distinct edges (distinct protein interactions present)
* Distinct Nodes (distinct proteins/drug present)
* Distribution of Edge scores in general and for each specific drug/protein
* Distribution of number of interaction per protein/drug 
* Amongst the shared proteins between the two networks, the difference and similarities in set of protein interactions present 

Because the scores in the STRING database were out of 1000,  ended up normalizing the STRING edge scores to a value between [0,1]. I removed any self edges and took the average scores for the redundant protein IDs. There are 4 important types of output folders saved when using the jupyter notebooks: NOT_FILTERED, FILTERED_FOR_SAME_EGDES, FILTERED_FOR_SAME_NODES, and FILTERED_FOR_CLOSE_NODES. For almost every ipython notebook, there is an output folder.


| Filter Type | Description |
| :---        |     :---    |
| NOT_FILTERED | No records/rows of protein-protein interactions were filtered out of the DataFrames for STRING or PathFX     |
| FILTERED_FOR_SAME_NODES | Only Proteins/Drugs (nodes in network graph) that were present in both STRING and PathFX were selected for analysis |
| FILTERED_FOR_SAME_EDGES | Only interactions (edges in network graph) that were present in both STRING and PathFX were selected for analysis |
| FILTERED_FOR_CLOSE_NODES | All the shared common proteins between PathFX and STRING + whatever nodes those common shared proteins connect to; PathFX and STRING will show a difference in smaller set of proteins - this as a "cleaner cut version" of the NOT_FILTER since STRING has so many other proteins not even present PathFX | 

The `perform_network_analysis.ipynb` executes all the notebooks and outputs everything in one execution. It takes 50 min for everything to be finished. It is also possible to run each filter analysis one-by-one instead if needed. 



### Built With

* [Anaconda](https://www.anaconda.com/)
* [NetworkX](https://networkx.org/)
* [Pandas](https://pandas.pydata.org/pandas-docs/stable/index.html)
* [Numpy](https://numpy.org/doc/stable/index.html)
* [Seaborn](https://seaborn.pydata.org/#)


<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Prerequisites

Technical Requirments:
* Must have at least 30GB of free memory on computer

List things you need to use the software and how to install them.
1. Install latest version of [Python](https://www.python.org/)
2. Install [Anaconda](https://www.anaconda.com/)
3. Install all the necessary packages, including [NetworkX](https://networkx.org/)
   ```sh
   pip install networkx
   ```

### Installation

1. Clone the repo
2. Download the entire inputs/ folder from the [network_analysis Google Drive Folder](https://drive.google.com/drive/u/1/folders/15Y8fLZutM89gyrMP6EVtuR6t_9Bc9HSi)
3. Open Anaconda environment and run JuypterLab
4. Open `perform_network_analysis.ipynb` on the JuypterLab Notebook


<!-- USAGE EXAMPLES -->
## Usage

Run all cells in the `perform_network_analysis.ipynb` on the JuypterLab Notebook and wait for execution to finish (Runtime ~50 min)

_For more information, please refer to the [Documentation](https://drive.google.com/drive/u/2/folders/1KlysvrHgyPZI8tK6X3x-Cjain11xNROu)_

<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/github_username/repo_name/issues) for a list of proposed features (and known issues).

<!-- LICENSE -->
## License

Distributed under the License. See `LICENSE` for more information.

<!-- CONTACT -->
## Contact

Arya Bakhtiar 
email: bakhtiararya [at] berkeley.edu

<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [GitHub Emoji Cheat Sheet](https://www.webpagefx.com/tools/emoji-cheat-sheet)
* [Img Shields](https://shields.io)
* [Choose an Open Source License](https://choosealicense.com)
* [GitHub Pages](https://pages.github.com)
* [Animate.css](https://daneden.github.io/animate.css)
* [Loaders.css](https://connoratherton.com/loaders)
* [Slick Carousel](https://kenwheeler.github.io/slick)
* [Smooth Scroll](https://github.com/cferdinandi/smooth-scroll)
* [Sticky Kit](http://leafo.net/sticky-kit)
* [JVectorMap](http://jvectormap.com)
* [Font Awesome](https://fontawesome.com)



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->

[contributors-shield]:  https://img.shields.io/github/contributors/aryastark5/network_analysis.svg?style=for-the-badge
[contributors-url]: https://github.com/aryastark5/network_analysis/graphs/contributors

[forks-shield]: https://img.shields.io/github/forks/aryastark5/network_analysis.svg?style=for-the-badge
[forks-url]: https://github.com/aryastark5/network_analysis/network/members

[stars-shield]: https://img.shields.io/github/stars/aryastark5/network_analysis.svg?style=for-the-badge
[stars-url]: https://github.com/aryastark5/network_analysis/stargazers

[issues-shield]: https://img.shields.io/github/issues/aryastark5/network_analysis.svg?style=for-the-badge
[issues-url]: https://github.com/aryastark5/network_analysis/issues

[license-shield]: https://img.shields.io/github/license/aryastark5/network_analysis.svg?style=for-the-badge
[license-url]: https://github.com/aryastark5/network_analysis/blob/main/LICENSE.txt

[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/
