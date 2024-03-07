# tecosa_demo_risk_accessment

## Generate a .tif file of a georeferenced map of KTH
Follow this [video](https://www.youtube.com/watch?v=RjkZgE_WVBk)

The generated .tif file can be found at KTH_geofig.tif

## Download road features from OpenStreetMap as .geojson using an API.
     curl --get 'https://osm.buntinglabs.com/v1/osm/extract' \
          --data "tags=highway=*" \
          --data "api_key=1kduLalBD857u3" \
          --data "bbox=18.068134,59.347679,18.071563,59.351561" \
          -o KTH_path_test.geojson

<ul>
  <li>enter your own api key</li>
  <li>defined bounding box (Lng / Lat)</li>
</ul>

<ul>
  <li>This step needs an API from: <a href="https://buntinglabs.com/solutions/openstreetmap-extracts">Bunting Labs</a></li>
  <li>Bounding box can be found at: <a href="http://bboxfinder.com/">bboxfinder</a></li>
  <li>Tags can be found at: <a href="https://taginfo.openstreetmap.org/">openstreetmap_tag</a></li>
</ul>

The generated .geojson file can be found at KTH_allroad.geojson
