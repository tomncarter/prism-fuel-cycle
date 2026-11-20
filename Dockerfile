FROM openmc/openmc:v0.15.3
RUN pip install --no-cache-dir pyarrow jupyterlab matplotlib pandas
