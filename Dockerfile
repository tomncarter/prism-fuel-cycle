FROM openmc/openmc:latest
RUN pip install --no-cache-dir pyarrow jupyterlab matplotlib pandas
