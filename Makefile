html:
	jupyter-book build .

conf.py: _config.yml _toc.yml
	jupyter-book config sphinx .

env:
	chmod a+x envscript.sh
	bash -ic envscript.sh
#mamba create -f environment.yml -p $HOME/envs/ligo
#conda activate ligo
#pip install pytest
#python2 -m ipykernel install --user --name ligo --display_name "Python2"
    

html-hub: conf.py
	sphinx-build  . _build/html -D html_baseurl=${JUPYTERHUB_SERVICE_PREFIX}/proxy/absolute/8000
	@echo "Start the Python http server and visit:"
	@echo "https://stat159.datahub.berkeley.edu/user-redirect/proxy/8000/index.html"

.PHONY: clean
clean:
	rm -rf _build/html/
	rm -rf figures
	rm -rf audio