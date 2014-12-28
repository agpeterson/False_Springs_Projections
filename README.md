# False Spring Projections

Spatiotemporal patterns of green-up dates, last spring freezes, and false springs are modeled across the contiguous United States from 1950 to 2099 using the MACAv2-METDATA dataset.

## Methodology
- For each model:
	- Calculate historical (1950-2005) mean value for each pixel.
	- Take difference between all annual values and historical mean.
	- Average over climatologies (2020-2049, others?).
- 