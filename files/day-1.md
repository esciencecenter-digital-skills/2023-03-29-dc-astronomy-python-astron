![](https://esciencecenter-digital-skills.github.io/2023-03-29-dc-astronomy-python-astron/assets/img/logos.png)


# Collaborative Document

2023-03-29 Astronomical Data Carpentry

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

### Day 1
| Time  | Topic                                                         |
| ----- | ------------------------------------------------------------- |
| 10:00 | Welcome / Coffee                                              |
| 10:15 | Intro / Basic Queries                                         |
| 11:45 | Coordinate Transformations (break 11:30 – 11:45)              |
| 13:00 | Lunch                                                         |
| 14:00 | Plotting and Tabular Data (break 15:00 - 15:15)               |
| 15:15 | Plotting and Pandas (break 16:00 – 16:15)                     |
| 16:30 | Buffer and flash-forward to tomorrow, feedback (borrel 17:30) |

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
trees

## 🎓 Certificate of attendance
If you attend the full workshop you can request a certificate of attendance by emailing to training@esciencecenter.nl .

## 🔧 Exercises
#### Excercise 1
One of the other tables we will use is `gaiadr2.panstarrs1_original_valid`. Use `load_table` to get the metadata for this table. How many columns are there and what are their names?

#### Solution 1
```python
# load the metadata
panstarrs_metadata = Gaia.load_table('gaiadr2.panstarrs1_original_valid')
print(panstarrs_metadata)

# print the column names
for column in panstarrs_metadata.columns:
    print(column.name)
```

#### Excerise 2

Read the [documentation](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html) of this table and choose a column that looks interesting to you. Add the column name to the query and run it again. What are the units of the column you selected? What is its data type?

#### Solution 2
``` python
query1_with_phot_g_mean_mag = """SELECT 
TOP 10
source_id, ra, dec, parallax, phot_g_mean_mag
FROM gaiadr2.gaia_source
"""
```

#### Exercise 3
The clauses in a query have to be in the right order. Go back and change the order of the clauses in ```query2``` and run it again. The modified query should fail, but notice that you don’t get much useful debugging information.

For this reason, developing and debugging ADQL queries can be really hard. A few suggestions that might help:
- Whenever possible, start with a working query, either an example you find online or a query you have used in the past.
- Make small changes and test each change before you continue.
- While you are debugging, use TOP to limit the number of rows in the result. That will make each test run faster, which reduces your development time.
- Launching test queries synchronously might make them start faster, too.

#### Excersise 4
[Read about SQL operators here](https://www.w3schools.com/sql/sql_operators.asp) and then modify the previous query to select rows where ```bp_rp``` is between ```-0.75``` and ```2```.

Previous query:
```python
query2 = """SELECT 
TOP 3000
source_id, ra, dec, pmra, pmdec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1
"""
```

#### Solution 4
```python
query2_sol1 = """SELECT 
TOP 3000
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1 
  AND bp_rp > -0.75 AND bp_rp < 2
"""
```
include bp_rp in the column selection to check yourself the if the where clause works as expected.
You can also use the BETWEEN operator instead; ``` AND bp_rp BETWEEN -0.75 AND 2 ```

#### Exercise 5
Create a quantity that represents 5 arcminutes and assign it to a variable called radius.

Then convert it to degrees.

#### Solution 5
```python
radius = 5 * u.arcmin
print(radius)

radius.to(u.degree)
```

#### Exercise 6
When you are debugging queries like this, you can use ```TOP``` to limit the size of the results, but then you still don’t know how big the results will be.

An alternative is to use ```COUNT```, which asks for the number of rows that would be selected, but it does not return them.

In the previous query, replace ```TOP 10 source_id``` with ```COUNT(source_id)``` and run the query again. How many stars has Gaia identified in the cone we searched?

#### Solution 6
``` python
count_cone_query = """SELECT 
COUNT(source_id)
FROM gaiadr2.gaia_source
WHERE 1=CONTAINS(
  POINT(ra, dec),
  CIRCLE(88.8, 7.4, 0.08333333))
"""

count_cone_job = Gaia.launch_job(count_cone_query)
count_cone_results = count_cone_job.get_results()
count_cone_results
```

#### Exercise 7
Find the location of GD-1 in ICRS coordinates.

- Create a ```SkyCoord``` object at 0°, 0° in the GD-1 frame.
- Transform it to the ICRS frame.

Hint: Because ICRS is a standard frame, it is built into Astropy. You can specify it by name, ```icrs``` (as we did with ```galactic```).

#### Solution 7
``` python
# create SkyCoord object
phi1 = 0.0 * u.degree
phi2 = 0.0 * u.degree
origin_gd1 = SkyCoord(phi1=phi1, phi2=phi2, frame=gd1_frame)

# transform to ICRS frame
origin_gd1.transform_to('icrs')
```

 #### <u> Day 1 afternoon </u>

#### Exercise 8
In the call to ```plt.plot```, use the keyword argument ```markersize``` to make the markers smaller.

Then add the keyword argument alpha to make the markers partly transparent.

Adjust these arguments until you think the figure shows the data most clearly.

Note: Once you have made these changes, you might notice that the figure shows stripes with lower density of stars. These stripes are caused by the way Gaia scans the sky, which [you can read about here](https://www.cosmos.esa.int/web/gaia/scanning-law). The dataset we are using, [Gaia Data Release 2](https://www.cosmos.esa.int/web/gaia/dr2), covers 22 months of observations; during this time, some parts of the sky were scanned more than others.

#### Solution 8
You can use for example:
```python
x = polygon_results['ra']
y = polygon_results['dec']
plt.plot(x, y, 'ko', markersize=0.1, alpha=0.1)

plt.xlabel('ra (degree ICRS)')
plt.ylabel('dec (degree ICRS)')
```

## 🧠 Collaborative Notes
### Day 1 morning

```python
from astroquery.gaia import Gaia
```

```python
# load all tables
tables = Gaia.load_tables(only_names=True)
```

```python
# print all tables that are available in the Gaia database
[t.name for t in tables]
```
#### Tables that will be used:
- gaiadr2.gaia_source
- gaiadr2.panstarrs1_original_valid
- gaiadr2.panstarrs1_best_neighbour

```python
table_metadata = Gaia.load_table('gaia_source')
```

```python
print(table_metadata)
```

```python
# look at the availbale columns in the data
for column in table_metadat.columns:
    print(column.name)
```

#### Writing queries
```python
# write an ASQL query
query1 = """SELECT
TOP 10
source_id, ra, dec, parallax
FROM gaiadr2.gaia_source
"""
```

```python
# send query to server
job1 = Gaia.launch_job(query1)
job1
```

```python
# show metadata of forthcoming results
print(job1)
```

```python
# get the results of the job
results1 = job1.get_results()
type(results1)
```

```python
# show the results in table form
results1
```

#### Asynchronous queries

```launch_job``` asks the server to run the job “synchronously”, which normally means it runs immediately. But synchronous jobs are limited to 2000 rows. For queries that return more rows, you should run “asynchronously”, using ```launch_job_async```, which mean they might take longer to get started.

``` python
# write a second query using WHERE
query2 = """SELECT 
TOP 3000
source_id, ra, dec, pmra, pmdec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1
"""
```
A WHERE clause indicates which rows we want; in this case, the query selects only rows “where” parallax is less than 1

``` python
# launch the query (asynchronously)
job2 = Gaia.launch_job_async(query2)
job2
```

``` python
# collect the results
results2 = job2.get_results()
results2
```

#### Formatting queries
``` python
query3_base = """SELECT
TOP 10
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2
"""
```

``` python
# define the columns to select
columns = ["source_id", "ra", "dec", "pmra", "pmdec", "parallax"]
```

``` python
query3 = query3_base.format(columns=", ".join(columns))
```

#### Coordinate Transformations

``` python
import astropy.units as u
```

``` python
# show all available units included in units
dir(u)
```

``` python
angle = 10 * u.degree
type(angle)
angle
```

``` python
# convert units
angle.to(u.arcmin)
```

If you add quantities, Astropy converts them to compatible units, if possible:
```python
angle + 30 * u.arcmin
```

```python
# if units are not compatible, you get an error
angle + 5 * u.kg
```


```python
# select a region
cone_query = """SELECT
TOP 10
source_id
FROM qaiadr2.gaia_source
WHERE CONTAINS(
    POINT(ra, dec),
    CIRCLE(88.8, 7.4, 0.83333333))"""
```
``` CONTAINS, POINT, CIRCLE``` are ADQL specific functions

(if you started a new notebook, import ```astroquery.gaia``` again)
``` python
from astroquery.gaia import Gaia

cone_job = Gaia.launch_job(cone_query)
cone_job
```

```python
cone_results = cone_job.get_results()
```

```python
from astropy.coordinates import SkyCoord

ra = 88.8 * u.degree
dec = 7.4 * u.degree
coord_icrs = SkyCoord(ra=ra, dec=dec, frame='icrs')
coord_icrs
```

```python
# convert the coordinates
coord_icrs.transform_to('galactic')
coord_icrs.ra
```

``` python
from gala.coordinates import GD1Koposov10

# define new coordinate frame
gd1_frame = GD1Koposov10()
gd1_frame
```

```python
# transform coordinates to gd1_frame
coord_icrs.transform_to(gd1_frame)
```

Select a rectangle
```python
phi1_min = -55 * u.degree
phi1_max = -45 * u.degree
phi2_min = -8 * u.degree
phi2_max = 4 * u.degree
```

```python
# write a function to create rectangles
def make_rectangle(x1, x2,y1, y2):
    xs = [x1, x1, x2, x2]
    ys = [y1, y2, y2, y1]
    return xs, ys
```

```python
phi1_rect, phi2_rect = make_rectangle(
    phi1_min, phi1_max, phi2_min, phi2_max)
```

```python
corners = SkyCoord(phi1=phi1_rect, phi2=phi2_rect, frame=gd1_frame)
corners
```

```python
list(zip(phi1_rect, phi2_rect))
```

```python
# transform to icrs coordinates
corners_icrs = corners.transform_to('icrs')
```

In order to use this polygon as part of an ADQL query, we have to convert it to a string with a comma-separated list of coordinates
```python
"""
POLYGON(143.65, 20.98, 
        134.46, 26.39, 
        140.58, 34.85, 
        150.16, 29.01)
"""
```

```python
# use to_string to produce a list of strings
corners_list_str = corners_icrs.to_string()
corners_list_str
```

```python
# to get the coordinates into ADQL format
corners_single_str = " ".join(corners_list_str)
corners_single_str.replace(' ', ', ')
```

```python
def skycoord_to_string(skycoord):
    """Convert a one-dimenstional list of SkyCoord to string for Gaia's query format."""
    corners_list_str = skycoord.transform_to('icrs').to_string()
    corners_single_str = ' '.join(corners_list_str)
    return corners_single_str.replace(' ', ', ')
```

```python
sky_point_list = skycoord_to_string(corners_icrs)
sky_point_list
```

```python
columns = ["source_id", "ra", "dec", "pmra", "pmdec", "parallax"]
```

```python
query3_base = """SELECT
TOP 10 
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
    AND bp_rp BETWEEN -0.75 AND 2
    AND 1 = CONTAINS(POINT(ra, dec), POLYGON({sky_point_list}))"""
```

```python
query3 = query3_base.format(columns=columns, sky_point_list=skycoord_to_string(corners_icrs)
print(query3)
```

```python
job3 = Gaia.launch_job_async(query3)
print(job3)
```

```python
# obtain and check length of results
results3 = job3.get_results()
len(results3)
```

```python
#save results to file
filename = 'gd1_results.fits'
results3.write(filename, overwrite=True)
```

### Day 1 afternoon
```python
import astropy.units as u
from astropy.coordinates import SkyCoord
from gaia.coordinates import GD1Koposov10
from astropy.table import Table

from episode_functions import *
```

```python
# load the saved results into a table
filename = "gd1_results.fits"
polygon_results = Table.read(filename)
```

```python
# check the data
len(polygon_results)
# or
polygon_results.info()
```

```python
# check out the column names
polygon_results.colnames
```

```python
# Select only the ra column
polygon_results['ra']
# Or specific rows from that column
polygon_results['ra'][0:3]
```

```python
# check the data type in the table
type(polygon_results["ra"][0])
```

```python
import matplotlib.pyplot as plt
```

```python
# to make sure your plots show in this notebook
%matpltlib inline
```

```python
x = polygon_results["ra"]
y = polygon_results["dec"]

plt.plot(x, y, "ro")
plt.xlabel("ra (degree ICRS)")
plt.ylabel("dec (degree ICRS)")
```

```python
gd1_frame = GD1Koposov10()
```

```python
skycoord = SkyCoord(ra=polygon_results["ra"], dec=polygon_results["dec"])
```

```python
distance = 8 * u.kpc
radial_velocity = 0 * u.km/y.s

skycoord = SkyCoord(ra=polygon_results['ra'], 
                    dec=polygon_results['dec'],
                    pm_ra_cosdec=polygon_results['pmra'],
                    pm_dec=polygon_results['pmdec'], 
                    distance=distance, 
                    radial_velocity=radial_velocity)
```

```python
# transform the coordinates to the gd1 frame
transformed = skycoord.transform_to(gd1_frame)
```

```python
from gala.coordinates import reflex_correct
```

```python
skycoord_gd1 = reflex_correct(transformed)
```

```python
# plot the transformed positions (in the gd1 frame)
x = skycoord_gd1.phi1
y = skycoord_gd1.phi2

plt.plot(x, y, 'ko', markersize=0.1, alpha=0.1)

plt.xlabel('phi1 (degree GD1)')
plt.ylabel('phi2 (degree GD1)')
```

```python
type(polygon_results)
type(skycoord_gd1)
```

```python
# add two columns to polygon results
polygon_results["phi1"] = skycoord_gd1.phi1
polygon_results["phi2"] = skycoord_gd1.phi2
# check out the new columns
polygon_results.info()
```

```python
# add the proper motions in phi1 and phi2
polygon_results["pm_phi1"] = skycoord_gd1.pm_phi1_cosphi2
polygon_results['pm_phi2'] = skycoord_gd1.pm_phi2
# check out the new columns
polygon_results.info()
```

#### Pandas dataframes
```python
import pandas

# load the results into a pandas dataframe
results_df = polygon_results.to_pandas()
```

```python
results_df.shape
```

```python
# show the top 5 rows of the dataframe
results_df.head()
```

```python
def make_dataframe(table):
    """Transform coordinates from ICRS to GD-1 frame.
    
    table: Astropy Table
    
    returns: Pandas DataFrame"""
    distance = 8 * u.kpc
    radial_velocity = 0 * u.km/u.s
    
    # Create a SkyCoord object with the coordinates and proper motions
    # in the input table
    skycoord = SkyCoord(
        ra=table['ra'], 
        dec=table['dec'],
        pm_ra_cosdec=table['pmra'],
        pm_dec=table['pmdec'],
        distance=distance,
        radial_velocity=radial_velocity
    )
    
    # Define the GD-1 reference frame
    gd1_frame = GD1Koposov10()

    # Transform input coordinates to the GD-1 reference frame
    transformed = skycoord.transform_to(gd1_frame)

    # Correct GD-1 coordinates for solar system motion around galactic center
    skycoord_gd1 = reflex_correct(transformed)

    #Add GD-1 reference frame columns for coordinates and proper motions
    table['phi1'] = skycoord_gd1.phi1
    table['phi2'] = skycoord_gd1.phi2
    table['pm_phi1'] = skycoord_gd1.pm_phi1_cosphi2
    table['pm_phi2'] = skycoord_gd1.pm_phi2

    # Create DataFrame
    df = table.to_pandas()
    
    return df
```

```python
results_df = make_dataframe(polygon_results)
```

```python
# Save the dataframe
filename = "gd1_data.hdf"
results_df.to_hdf(filename, "results_df", mode="w")
```
the dataframe will be called "results_df" within the hdf file.


```python
# Plot the proper motion
x = results_df['pm_phi1']
y = results_df['pm_phi2']

plt.plot(x, y, 'ko', markersize=0.3, alpha=0.3)

plt.xlabel('Proper motion phi1')
plt.ylabel('Proper motion phi2')

# zoom in
plt.xlim(-12, 8)
plt.ylim(-10, 10)
```

```python
phi2_min = -1.0 * u.degree
phi2_max = 1.0 * u.degree

phi2 = results_df["phi2"]
type(phi2)
```

```python
# create a mask to make a subselection of the stars
mask = (phi2 > phi2_min)
type(mask)
```

```python
mask = (phi2 > phi2_min) & (phi2 < phi2_max)
```

```python
mask.sum()
```

```python
# apply the mask to the dataframe
centerline_df = results_df[mask]
```

```python
# Define a function to plot the proper motion
def plot_proper_motion(df):
    """Plot proper motion for input dataframe"""
    x = df['pm_phi1']
    y = df['pm_phi2']

    plt.plot(x, y, 'ko', markersize=0.3, alpha=0.3)

    plt.xlabel('Proper motion phi1')
    plt.ylabel('Proper motion phi2')

    # zoom in
    plt.xlim(-12, 8)
    plt.ylim(-10, 10)
```

``` python
plot_proper_motion(centerline_df)
```

``` python
pm1_min = -8.9
pm1_max = -6.9
pm2_min = -2.2
pm2_max = 1.0
```

``` python
pm1_rect, pm2_rect = make_rectangle(pm1_min, pm1_max, pm2_min, pm2_max)
```

``` python
plt_proper_motion(centerline_df)
# overplot created rectangle
plt.plot(pm1_rect, pm2_rect, '-')
```

``` python
def between(series, low, high):
    return (series > low) & (series < high)
```

``` python
pm1 = results_df["pm_phi1"]
pm2 = results_df["pm_phi2"]
pm_mask = (between(pm1, pm1_min, pm1_max) & between(pm2, pm2_min, pm2_max))
```

``` python
pm_mask.sum()
```

``` python
# apply created mask to results_df
selected_df = results_df[pm_mask]
```

``` python
x = selected_df['phi1']
y = selected_df['phi2']
plt.plot(x, y, 'ko', markersize=0.3, alpha=0.3)

plt.xlabel('phi1 (degree GD1)')
plt.ylabel('phi2 (degree GD1)')
```

``` python
def plot_pm_selection(df)):
    x = df['phi1']
    y = df['phi2']

    plt.plot(x, y, "ko", ms=0.3, alpha=0.3)

    plt.xlabel('$\phi_1$ (deg)')
    plt.ylabel('$\phi_2$ (deg)')
    
    plt.title('proper motion selection', fontsize='medium')
    plt.axis("equal")
```

``` python
# save the dataframe
selected_df.to_hdf(filename, 'selected_df')
```

```python
# save centerline_df
centerline_df.to_hdf(filename, "centerline_df")
```

```python
# show what is in the file
with df.HDFStore(filename) as hdf:
    print(hdf.keys())
```

## 📚 Resources
[Paper](https://iopscience.iop.org/article/10.3847/2041-8213/aad7b5) by Price-Whelan and Bonaca
[astroquery](https://astroquery.readthedocs.io/en/latest/)
[astroquery.gaia](https://astroquery.readthedocs.io/en/latest/gaia/gaia.html)
[Gaia documentation](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html) on the table gaia_source
[ADQL documentation](https://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.pdf)

## Feedback

### To improve

### Like
