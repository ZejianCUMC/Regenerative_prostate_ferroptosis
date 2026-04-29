import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score



def predict_and_save_results(model_path, test_data_path, output_path):
    # Load the trained model
    testmodel = xgb.XGBClassifier()
    testmodel.load_model(model_path)

    # Load test data
    test_data = pd.read_csv(test_data_path)

    # Separate true labels from features
    true_labels = test_data['CellType']
    test_data = test_data.drop('CellType', axis=1)

    # Encode the true labels
    label_encoder = LabelEncoder()
    label_encoder.fit(true_labels)

    # Handle missing columns
    missing_cols = set(testmodel.feature_names_in_) - set(test_data.columns)
    missing_cols_df = pd.DataFrame(0, index=test_data.index, columns=list(missing_cols))
    test_data = pd.concat([test_data, missing_cols_df], axis=1)

    # Reorder columns to match the training data
    X_combined = pd.DataFrame(columns=testmodel.feature_names_in_)
    test_data = test_data[X_combined.columns]

    # Predict labels
    predicted_labels = testmodel.predict(test_data)
    predicted_labels_transformed = label_encoder.inverse_transform(predicted_labels)

    # Predict probabilities
    probs = testmodel.predict_proba(test_data)
    prob_df = pd.DataFrame(probs, columns=label_encoder.classes_)
    prob_df['Predicted_Class'] = predicted_labels_transformed
    prob_df['True_Class'] = true_labels

    # Save the results to a CSV file
    prob_df.to_csv(output_path, index=False)

# Example usage:
model_path = 'Intact.xgboost_model.json'
test_data_path = 'MW3_luminal_scaled.csv'
output_path = 'MW3_luminal_prediction.csv'

predict_and_save_results(model_path, test_data_path, output_path)
