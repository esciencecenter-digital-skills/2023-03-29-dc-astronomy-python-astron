![](https://i.imgur.com/iywjz8s.png)


# Collaborative Document

2023-03-30 Astronomical Data Carpentry

Welcome to The Workshop Collaborative Document.

This Document is synchronized as you type, so that everyone viewing this page sees the same text. This allows you to collaborate seamlessly on documents.

----------------------------------------------------------------------------

This is the Document for today: [link](<url>)

Collaborative Document day 1: [link](https://codimd.carpentries.org/RP0Ns-wVR9eLjQzicpUoVw#)

Collaborative Document day 2: [link](https://codimd.carpentries.org/E0BWinWsTpGEE---o08BVA#)

## 👮Code of Conduct

Participants are expected to follow these guidelines:
* Use welcoming and inclusive language.
* Be respectful of different viewpoints and experiences.
* Gracefully accept constructive criticism.
* Focus on what is best for the community.
* Show courtesy and respect towards other community members.
 
## ⚖️ License

All content is publicly available under the Creative Commons Attribution License: [creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/).

## 🙋Getting help

To ask a question, just raise your hand.

If you need help from a helper, place a pink post-it note on your laptop lid. A helper will come to assist you as soon as possible.

## 🖥 Workshop website

- [Workshop website](https://esciencecenter-digital-skills.github.io/2023-03-29-dc-astronomy-python-astron/)
- [🛠 Setup instructions](https://datacarpentry.org/astronomy-python/setup.html)
- [Download files](https://figshare.com/ndownloader/files/35777540)

## 👩‍🏫👩‍💻🎓 Instructors
Johan Hidding, Hanno Spreeuw

## 🧑‍🙋 Helpers
Laura Ootes

## 👩‍💻👩‍💼👨‍🔬🧑‍🔬🧑‍🚀🧙‍♂️🔧 Roll Call
Name/ pronouns (optional) / job, role / social media (twitter, github, ...) / background or interests (optional) / city

## 🗓️ Agenda

### Day 2
| Time  | Topic                                     |
| ----- | ----------------------------------------- |
| 9:00  | Recap, Transform and Select (break 10:30) |
| 10:45 | Join                                      |
| 12:30 | Lunch                                     |
| 13:30 | Photometry (break 14:30 – 14:45)          |
| 14:45 | Visualization (break 15:30 – 15:45)       |
| 16:30 | Recap, questions, feedback (finish 17:00) |


## 🏢 Location logistics
:deciduous_tree: 🪺 :evergreen_tree: :mushroom: 🍃 :evergreen_tree: 🪶 :seedling: :spider_web: :deciduous_tree: :fallen_leaf: :chestnut: 
:owl: :mouse: :fox_face: :rabbit: :deer: :eagle: :wolf: :ox: :boar: 🐑 :chipmunk: :hedgehog: :bat: 
:beetle: :ant: 🪱 :spider: :bee: :snail: :butterfly: 
## 🎓 Certificate of attendance
If you attend the full workshop you can request a certificate of attendance by emailing to training@esciencecenter.nl .

## 🔧 Exercises
### Day 2 morning
#### Exercise 1
Looking at the proper motion of the stars we identified along the centerline of GD-1, in the ICRS reference frame define a rectangle (```pmra_min```, ```pmra_max```, ```pmdec_min```, and ```pmdec_max```) that encompass the proper motion of the majority of the stars near the centerline of GD-1 without including to much contamination from other stars.

#### Solution 1
```python
pmra_max = selected_df.pmra.max()
pmra_min = selected_df.pmra.min()
pmdec_max = selected_df.pmdec.max()
pmdec_min = selected_df.pmdec.min()
pmra_rect, pmdec_rect = make_rectangle(pmra_min, pmra_max, pmdec_min, pmdec_max)
```

#### Exercise 2
Define ```candidate_coord_pm_query_base```, starting with ```candidate_coord_query_base``` and adding two new ```BETWEEN``` clauses to select stars whose coordinates of proper motion, pmra and pmdec, fall within the region defined by ```pmra_min```, ```pmra_max```, ```pmdec_min```, and ```pmdec_max```. In the next exercise we will use the format statement to fill in the values we defined above.

#### Solution 2
```python
candidate_coord_query_base = """SELECT
    {columns}
    FROM gaiadr2.gaia_source
    WHERE parallax < 1
        AND bp_rp BETWEEN -0.75 AND 2
        AND 1 = CONTAINS(POINT(ra,dec),
                         POLYGON({sky_point_list}))
        AND pmra BETWEEN {pmra_min} AND  {pmra_max}
        AND pmdec BETWEEN {pmdec_min} AND {pmdec_max}
"""
```


#### Exercise 3
Now we are ready to bring in the Pan-STARRS table. Starting with the previous query, add a second `JOIN` clause that joins with `gaiadr2.panstarrs1_original_valid`, gives it the abbreviated name ps, and matches `original_ext_source_id` from the best neighbor table with obj_id from the Pan-STARRS table.
Add `g_mean_psf_mag` and `i_mean_psf_mag` to the column list, and run the query. The result should contain 490 rows and 9 columns.

#### Solution 3

```python
# join the tables
cone_base_query = """
SELECT
{columns}
FROM gaiadr2.gaia_source AS gaia
JOIN gaiadr2.panstarrs1_best_neighbour AS best
    ON gaia.source_id = best.source_id
JOIN gaiadr2.panstarrs1_original_valid as ps
  ON best.original_ext_source_id = ps.obj_id
    -- USING(source_id)
WHERE 1 = CONTAINS(
    POINT(gaia.ra, gaia.dec)
    CIRCLE(88.8, 7.4, 0.0833333)
)"""
```


```python
columns = ["gaia.source_id", 
           "gaia.ra", 
           "gaia.dec", 
           "gaia.pmra", 
           "gaia.pmdec", 
           "best.best_neighbour_multiplicity", 
           "best.number_of_mates",
           "ps.g_mean_psf_mag",
           "ps.i_mean_psf_mag"
          ]
cone_query = cone_base_query.format(
    columns = ". ".join(columns))
print(cone_query)
```

#### Exercise 4
Create a new query base called `candidate_join_query_base` that combines the `WHERE` clauses from the previous query with the JOIN clauses for the best neighbor and Pan-STARRS tables. Format the query base using the column names in `column_list`, and call the result `candidate_join_query`.

Hint: Make sure you use qualified column names everywhere!

Run your query and download the results. The table you get should have 4300 rows and 9 columns.

### Day 2 afternoon
#### Exercise 5
When we encounter a new object, it is good to create a toy example to test that it does what we think it does. Define a list of two points (represented as two tuples), one that should be inside the polygon and one that should be outside the polygon. Call contains_points on the polygon we just created, passing it the list of points you defined, to verify that the results are as expected.

#### Solution 5
For example:
```python
polygon.contains_points((0.4, 20))
```

#### Exercise 6
Boolean values are stored as 0s and 1s. `FALSE = 0` and `TRUE = 1`. Use this information to determine the number of stars that fall inside the polygon.

#### Solution 6
```python
inside_mask = polygon.contains_points(cmd_df)
inside_mask.sum()
# should return 640
#fraction:
candidate_df.shape[0] / inside_mask.sum()
```

#### Exercise 7
Think about the following questions:

1. What is the primary scientific result of this work?

2. What story is this figure telling?
3. In the design of this figure, can you identify 1 or 2 choices the authors made that you think are effective? Think about big-picture elements, like the number of panels and how they are arranged, as well as details like the choice of typeface.
4. Can you identify 1 or 2 elements that could be improved, or that you might have done differently?

#### Exercise 8
Plot the selected stars in `winner_df` using the `plot_cmd_selection` function and then choose any or all of these features and add them to the figure:

- To draw vertical lines, see [`plt.vlines`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.vlines.html) and [`plt.axvline`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axvline.html).
- To add text, see [`plt.text`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html).
- To add an annotation with text and an arrow, see [`plt.annotate`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.annotate.html).

Here is some [additional information about text and arrows](https://matplotlib.org/stable/tutorials/text/annotations.html).

## 🧠 Collaborative Notes
### Day 2 morning
```python
import astropy.units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from gala.coordinates import GD1Koposov10, GD1, reflex_correct
import matplotlip.pyplot as plt
import pandas as pd

from episode_functions import *
```
---  intermezzo ---
---
```python
filename = "gd1_data.hdf"
centerline_df = pd.read_hdf(filename, "centerline_df")
selected_df = pd.read_hdf(filename, "selected_df")
```

```python
#class rectangle
#    def __init__(self, pm1_min, pm1_max, pm2_min, pm2_max)
#        self.pm1...

@dataclass
class Rectangle:
    x_min: Quantity
    x_max: Quantity
    y_min: Quantity
    y_max: Quantity
        
    def sky_coords(self):
        xs = [self.x_min, self.x_min, self.x_max, self.x_max]
        ys = [self.y_min, self.y_max, self.y_max, self.y_min]
        return SkyCoord(xs, ys, frame=frame)
```

```python
rect = Rectangle(-8.9 * u.degee, -6.9 * u.degee, -2.2 * u.degee, 1.0 * u.degee)
```

```python
pm_phi1, pm_phi2 = rect.rect()

rect.sky_coords(gd1_frame).transform_to("icrs")
```

--- end of intermezzo ---
---

```python
pm1_min = -8.9
pm1_max = -6.9
pm2_min = -2.2
pm2_max = 1.0

pm1_rect, pm2_rect = make_rectangle(
    pm1_min, pm1_max, pm2_min, pm2_max)

gd1_frame = GD1Koposov10()
```

```python
# using the plot_proper_motion function defined yesterday
# Plot proper motion for the stars we selected along the center line of GD-1
# The selected rectangle
# the stars inside the rectange highlighted in green
plot_proper_motion(centerline_df)
plt.plot(pm1_rect, pm2_rect)

x = selected_df['pm_phi1']
y = selected_df['pm_phi2']

plt.plot(x, y, 'gx', markersize=0.3, alpha=0.3);
```

```python
pmra_max = selected_df.pmra.max()
pmra_min = selected_df.pmra.min()
pmdec_max = selected_df.pmdec.max()
pmdec_min = selected_df.pmdec.min()
pmra_rect, pmdec_rect = make_rectangle(pmra_min, pmra_max, pmdec_min, pmdec_max)
```

```python
# make the same plot in the icrs frame, which are stored in the columns pmra and pmdec
x = centerline_df.pmra
y = centerline_df.pmdec
plt.plot(x, y, 'ko', markersize=0.3, alpha=0.3)

x = selected_df.pmra
y = selected_df.pmdec
plt.plot(x, y, 'gx', markersize=1, alpha=0.3)
# the more specific rectangle
plt.plot(pmra_rect, pmdec_rect)

plt.xlabel('Proper motion ra (ICRS frame)')
plt.ylabel('Proper motion dec (ICRS frame)')

plt.xlim([-10, 5])
plt.ylim([-20, 5]);
```

Query the database
```python
candidate_coord_query_base = """SELECT {columns}
    FROM gaiadr2.gaia_source
    WHERE parallax < 1
        AND bp_rp BETWEEN -0.75 AND 2
        AND 1 = CONTAINS(POINT(ra,dec),
                         POLYGON({sky_point_list}))
"""
```

```python
# Define new rectangle
phi1_min = -70 * u.degree
phi1_max = -20 * u.degree
phi2_min = -5 * u.degree
phi2_max = 5 * u.degree

phi1_rect, phi2_rect = make_rectangle(
    phi1_min, phi1_max, phi2_min, phi2_max)
```

```python
corners = SkyCoord(phi1=phi_rect,
                   phi2=phi2_rect,
                   frame=gd1_frame)
corners_icrs = corners.transform_to("icrs")
```

```python
sky_point_list = skycoord_to_string(corners_icrs)
```

```python
# check your list
sky_point_list
```

```python
columns = ["source_id", "ra", "dec", "pmra", "pmdec"]
```

```python
candidate_coord_query = candidate_coord_query_base.format(
    columns=", ".join(columns),
    sky_point_list=sky_point_list)
print(candidate_coord_query)
```

```python
candidate_coord_query = candidate_coord_query_base.format(
    columns=", ".join(columns),
    sky_point_list=sky_point_list,
    pmra_min=pmra_min,
    pmra_max=pmra_max,
    pmdec_min=pmdec_min,
    pmdec_max=pmdec_max)
print(candidate_coord_pm_query)
```

```python
# launch job async
cc_job = Gaia.launch_job_async(candidate_coord_query)
print(cc_job)
```

```python
# get the resuls
cc_table = cc_job.get_results()
cc_results = cc_job.results()
cc_df = cc_results.to_pandas
```

```python
x = cc_df.ra
y = cc_df.dec
plt.plot(x, y, 'ko', markersize=0.3, alpha=0.3)
plt.xlabel("ra (ICRS)")
plt.ylabel("dec (ICRS)")
```

```python
# convert to gd1 frame
cc_gd1 = make_dataframe(cc_results)
```

```python
# plot proper motion selection using plot function defined yesterday
plt_pm_selection(cc_gd1)
```

### Join
```python
# load the data
ps_best_neighbour_meta = Gaia.load_table("gaiadr2.panstarrs1_best_neighbour")
```

```python
print(ps_best_neighbour_meta)
```

```python
# check out the available columns
for column in ps_best_neighbour_meta.columns:
    print(column.name)
```

```python
# set up the query
ps_best_neighbour_query = """ SELECT
TOP 5
source_id, best_neighbour_multiplicity, number_of_mates, original_ext_source_id
FROM gaiadr2.panstarrs1_best_neighbour
"""
```

```python
ps_best_neighbour_job = Gaia.launch_job(ps_best_neighbour_query)
ps_best_neighbout_job.results
```

```python
ps_valid_meta = Gaia.load_table("gaiadr2.panstarrs1_original_valid")
print(ps_valid_meta)
```

```python
for column in ps_valid_meta.columns:
    print(column.name)
```

```python
ps_valid_query = """SELECT
TOP 5
obj_id, g_mean_psf_mag, i_mean_psf_mag 
FROM gaiadr2.panstarrs1_original_valid"""

ps_valid_job = Gaia.launch_job(ps_valid_query)
```

```python
ps_valid_job.results
```

```python
# join the tables
cone_base_query = """
SELECT
-- TOP 10
{columns}
FROM gaiadr2.gaia_source as gaia

-- cross-table
JOIN gaiadr2.panstarrs1_best_neighbour as best
  ON gaia.source_id = best.source_id

-- pannstrarrs
JOIN gaiadr2.panstarrs1_original_valid as ps
  ON best.original_ext_source_id = ps.obj_id
WHERE 1=CONTAINS(
  POINT(gaia.ra, gaia.dec),
  CIRCLE(88.8, 7.4, 0.08333333))"""
```


```python
columns = ["gaia.source_id", 
           "gaia.ra", 
           "gaia.dec", 
           "gaia.pmra", 
           "gaia.pmdec", 
           "best.best_neighbour_multiplicity", 
           "best.number_of_mates",
           "ps.g_mean_psf_mag",
           "ps.i_mean_psf_mag"
          ]
cone_query = cone_base_query.format(
    columns = ". ".join(columns))

print(cone_query)
```

```python
cone_job = Gaia.launch_job_async(cone_query)
```

```python
cone_job.results
```

```python
candidate_coord_query_base = """SELECT {columns}
    FROM gaiadr2.gaia_source
    WHERE parallax < 1
        AND bp_rp BETWEEN -0.75 AND 2
        AND 1 = CONTAINS(POINT(ra,dec),
                         POLYGON({sky_point_list}))
    AND pmra BETWEEN {pmra_min} AND  {pmra_max}
    AND pmdec BETWEEN {pmdec_min} AND {pmdec_max}                  
```

```python
join_query_base = """
SELECT
-- TOP 10
{columns}
FROM gaiadr2.gaia_source as gaia

-- cross-table
JOIN gaiadr2.panstarrs1_best_neighbour as best
  ON gaia.source_id = best.source_id

-- pannstrarrs
JOIN gaiadr2.panstarrs1_original_valid as ps
  ON best.original_ext_source_id = ps.obj_id
 
-- deselect foreground
WHERE gaia.parallax < 1 
    AND gaia.bp_rp BETWEEN -0.75 AND 2
    -- areas on the sky
    AND 1=CONTAINS(
        POINT(gaia.ra, gaia.dec),
          CIRCLE(88.8, 7.4, 0.08333333))
    -- select on proper motion
    AND gaia.pmra BETWEEN {pmra_min} AND  {pmra_max}
    AND gaia.pmdec BETWEEN {pmdec_min} AND {pmdec_max}
  """
```

```python
columns = ['gaia.source_id',
            'gaia.ra',
            'gaia.dec',
            'gaia.pmra',
            'gaia.pmdec',
            'best.best_neighbour_multiplicity',
            'best.number_of_mates',
            'ps.g_mean_psf_mag',
            'ps.i_mean_psf_mag']

join_query = join_query_base.format(columns=', '.join(column))
```

```python
join_job = Gaia.launch_job_async(join_query)

```

```python
candidate_table = join_job.results
```

```python
# check
all(candidate_table["best_neighbour_multiplicity"]==1)
```

```python
sum(candidate_table["number_of_mates"] != 0)
```

```python
# other way to check the results
multiplicity = pd.Series(candidate_table['best_neighbour_multiplicity'])
multiplicity.describe()
```

```python
mates = pd.Series(candidate_table['number_of_mates'])
mates.describe()
```

```python
candidate_df = make_dataframe(candidate_table)
```

```python
# quick inspection of results
candidate_df.plot("phi1", "phi2", kind='hexbin', gridsize=50)
```

```python
# save the results
filename = "gd1_data.hdf"
candidate_df.to_hdf(filename, "candidate_df")
```
```python
from os.path import getsize

MB = 1024 * 1024
getsize(filename) / MB
```


### Day 2 afternoon

```python
# Plot colour-magnitude diagram
x = candidate_df['g_mean_psf_mag'] - candidate_df['i_mean_psf_mag']
y = candidate_df['g_mean_psf_mag']

plt.plot(x, y, 'ko', markersize=0.3, alpha=0.3)
plt.xlim([0, 1.5])
plt.ylim([22, 14])
plt.ylabel('Magnitude (g)')
plt.xlabel('Colour (g-i)')
```

```python
def plot_cmd(dataframe):
    x = dataframe['g_mean_psf_mag'] - dataframe['i_mean_psf_mag']
    y = dataframe['g_mean_psf_mag']

    plt.plot(x, y, 'ko', markersize=0.3, alpha=0.3)
    plt.xlim([0, 1.5])
    plt.ylim([22, 14])
    plt.ylabel('Magnitude (g)')
    plt.xlabel('Colour (g-i)')
```

```python
plot_cmd(candidate_df)
```

```python
filename = 'gd1_isochrone.hdf5'
iso_df = pd.read_hdf(filename, 'iso_df')
iso_df.head()
```

```python
plot_cmd(candidate_df)
plt.plot(iso_df['color_g_i'], iso_df['mag_g'])
```

```python
g_all = iso_df['mag_g']

g_mask = (g_all > 18.0) & (g_all < 21.5)
# check number of entries
g_mask.sum()
```

```python
# apply the mask
iso_masked = iso_df[g_mask]
iso_masked.head()
```

```python
g = iso_masked['mag_g']
left_color = iso_masked['color_g_i'] - 0.06
right_color = iso_masked['color_g_i'] + 0.12
```

```python
plot_cmd(candidate_df)

plt.plot(left_color, g, label='left colour')
plt.plot(right_color, g, label='right colour')

plt.legend()
```

```python
reverse_right_color = right_color[::-1]
```

```python
import numpy as np
color_loop = np.append(left_color, reverse_right_color)
color_loop.shape
```

```python
mag_loop = np.append(g, g[::-1])
mag_loop.shape
```

```python
plot_cmd(candidate_df)
plt.plot(color_loop, mag_loop)
```

```python
# make pandas dataframe out of loops
loop_df = pd.DataFrame()
loop_df['color_loop'] = color_loop
loop_df['mag_loop'] = mag_loop
loop_df.head()
```

```python
from matplotlib.patches import Polygon

polygon = Polygon(loop_df)
polygon
```

```python
filename = 'gd1_data.hdf'
loop_df.to_hdf(filename, 'loop_df')
```

```python
cmd_df = pd.DataFrame()
cmd_df['color'] = candidate_df['g_mean_psf_mag'] - candidate_df['i_mean_psf_mag']
cmd_df['mag'] = candidate_df['g_mean_psf_mag']

cmd_df.head()
```

```python
polygon.contains_points(cmd_df)
```

```python
winner_df = candidate_df[inside_mask]
```

```python
plot_cmd(candidate_df)
plt.plot(iso_df['color_g_i'], iso_df['mag_g'])
plt.plot(color_loop, mag_loop)

x = winner_df['g_mean_psf_mag'] - winner_df['i_mean_psf_mag']
y = winner_df['g_mean_psf_mag']

plt.plot(x, y, 'go', markeersize=0.5, alpha=0.5)
```

```python
def plot_cmd_selection(df):
    x = df['phi1']
    y = df['phi2']
    
    plt.plot(x, y, 'ko', markersize=0.7, alpha=0.9)
    
    plt.xlabel('$\phi_1$ [deg]')
    plt.ylabel('$\phi_2$ [deg]')
    plt.title('Proper motion + HR-diagram selection', fontsize='medium')
    plt.axis('equal')
```

```python
fig = plt.figure(figsize=(10,2.5))
plot_cmd_selection(winner_df)
```

```python
# add dataframe to hdf file
filename = 'gd1_data.hdf'
winner_df.to_hdf(filename, "winner_df")
```

```python
fig = plt.figure(figsize=(10,2.5))
ax = fig.add_subplot(1,1,1)
ax.tick_params(direction='in')
plot_cmd_selection(winner_df)

plt.axvline(-55, ls='--', color='gray', alpha=0.4, 
            dashes=(0.4), lw=2)
plt.text(-60, 5.5, "Previously \n  undetected", 
         fontsize="small", ha="right", va='top')

arrowprops = dict(color='gray', shrink=0.05, width=1.5,
                  headwidth=6, headlength=8, alpha=0.4)

plt.annotate('Spur', xy=(-33, 2), xytext=(-35, 5.5),
             arrowprops=arrowprops,
             fontsize='small')

plt.annotate('Gap', xy=(-22, -1), xytext=(-25, -5.5),
             arrowprops=arrowprops,
             fontsize='small');
```

```python
plt.rcParams['font.size'] = 14
```

```python
# check location of default parameter settings file
import matplotlib as mpl
mpl.matplotlib_fname()

```

```python
# Show available stylesheets
plt.style.available
```

```python
# Use specific stylesheet
plt.style.use('seaborn-paper')
```

```python
plt.style.use('./az-paper-twocol.mplstyle')

plot_cmd(candidate_df)

plt.plot(left_color, g, label='left color')
plt.plot(right_color, g, label='right color')

plt.legend();
```

## 📚 Resources
[Paper](https://iopscience.iop.org/article/10.3847/2041-8213/aad7b5) by Price-Whelan and Bonaca
[astroquery](https://astroquery.readthedocs.io/en/latest/)
[astroquery.gaia](https://astroquery.readthedocs.io/en/latest/gaia/gaia.html)
[Gaia documentation](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html) on the table gaia_source
[ADQL documentation](https://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.pdf)
[MIST](http://waps.cfa.harvard.edu/MIST/)
[AAS graphics guide](https://journals.aas.org/graphics-guide/)


## Feedback
- [Post-workshop survey](https://www.surveymonkey.com/r/82Y3JBT)
### To improve


### Like

